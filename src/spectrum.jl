import Statistics
using StaticArrays
using CPUTime
using Interpolations

using PhysConst
using SimpleLog

# Save convolved absorption and emission spectra to hdf5 file
function save_convolved_to_hdf5(λb, κb, ϵb, hdf5_path)
    groups = Dict("wl_bin" => Dict("wl" => λb, "kappa" => κb, "epsilon" => ϵb))
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end

# Save convolved absorption and emission spectra to ny file
function save_spectrum_as_npy(path, λl, κl, ϵl)
    groups = Dict("wl_lines" => Dict("wl" => λl, "kappa" => κl, "epsilon" => ϵl))
end

#function save_spectrum_as_hdf5(path, spec::Spectrum)
#    nb_λ = size(spec.λ, 1)
#    A = Matrix{Float64}(undef, 3, nb_λ)
#    for i in 1:nb_λ
#        A[1,i] = spec.λ[i];
#        A[2,i] = spec.κ_c[i];
#        A[3,i] = spec.ϵ_c[i];Planck_Ts
#    nb_zsteps = size(spec.z,1)
#    for istep in 1:nb_zsteps
#        p = y_at(spec.p_vs_h, spec.z[istep])
#        T = y_at(spec.t_vs_h, spec.z[istep])
#        N = p / (c_kB * T)
#        @info @sprintf("z = %8.2e, p = %8.2e, T = %8.2e, N = %8.2e", spec.z[istep], p, T, N)
#    end
#end

function compute_line_emission_and_absorption_iλ(ldspecies::LineData, Qinpspecies, cspecies, mspecies, T, N, p, T_ref, iλ)
    dΩ = 1.0
    β  = 1.0/(c_kB * T)

    #                  1      2   3     4  5     6    7    8    9    10   11   12    13    14    15     16
    #lines[i] = SA_F64[λ_ul0, E_l, E_u, S, A_ul, γ_a, γ_s, n_a, δ_a, g_u, g_l, B_ul, B_lu, ΔλL0, iso_m, iso_c]

    λ210  = ldspecies.λ210[iλ]
    E1    = ldspecies.E1[iλ]
    E2    = ldspecies.E2[iλ]
    A21   = ldspecies.A21[iλ]
    γair  = ldspecies.γair[iλ]
    γself = ldspecies.γself[iλ]
    nair  = ldspecies.nair[iλ]
    δair  = ldspecies.δair[iλ]
    g2    = ldspecies.g2[iλ]
    g1    = ldspecies.g1[iλ]
    B12   = ldspecies.B12[iλ]
    B21   = ldspecies.B12[iλ]
    mid   = ldspecies.mid[iλ]
    lid   = ldspecies.lid[iλ]+1

    Q     = Qinpspecies[lid](T)
    λ21   = λ210 / (1.0 + λ210 * δair * p)

    # γ = (Tref/T)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})
    dT = (T_ref/T)^nair
    γ = dT * (γair * p * (1.0 - cspecies) + γself * p * cspecies)

    ΔλL = λ21^2 * γ
    ΔλG = sqrt(2.0 * c_kB * T / mspecies[lid]) / c_c * λ21

    N1  = g1 * exp(- E1 * β) / Q * cspecies * N
    N2  = g2 * exp(- E2 * β) / Q * cspecies * N

    # ϵ_λ = h * c / λ_0 / (4 * π) * N_u * A_ul * f_λ
    ϵ = c_h * c_c / λ21 *  N2 * A21 * dΩ / (4.0 * π)

    # κ_λ = h / λ_0 * N_l * B_lu * (1 - N_u/N_l * g_l/g_u) * λ_0**2 / c * f_λ
    κ = c_h * λ21 / c_c * (N1 * B12 - N2 * B21)

    mid, lid, λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ
end

@doc raw"""
    Compute emission and absorption coefficients of the lines

    T - temperature
    N - atmosphere density
    p - pressure
    NCO2 - CO2 concentration

$γ = \left(\dfrac{T_{ref}}{T}\right)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})

Niso = N * NCO2 * iso_c

λ_ul = λ_ul0 / (1.0 + λ_ul0 * δ_a * p)

# γ = (Tref/T)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})
dT = (par.T_ref/T)^n_a
γ = dT * (γ_a * p * (1.0 - NCO2) + γ_s * p * NCO2)

ΔλL = λ_ul^2 * γ
ΔλG = sqrt(2.0 * c_kB * T / iso_m) / c_c * λ_ul
ΔλL_mean[iλ] = ΔλL
ΔλD_mean[iλ] = ΔλG

$N_l  = \dfrac{g_l}{Q(T, iso)} \exp(- E_l  β)  N_{iso}$
$N_u  = \dfrac{g_u}{Q(T, iso)} \exp(- E_u  β)  N_{iso}$

$ϵ(λ) = \dfrac{h c}{λ_0} N_u A_{ul} * f(λ) \dfrac{dΩ}{4 π}$
$κ(λ) = \dfrac{h c}{λ_0} N_l B_{lu}  \left(1 - \dfrac{N_u}{N_l}  \dfrac{g_l}{g_u}\right)  \dfrac{λ_0^2}{c} f(λ)$

compute_lines_emission_and_absorption(moleculardata, linedata, Nmolecules, T, N, p)
"""
function compute_lines_emission_and_absorption(moleculardata::Vector{MolecularData}, linedata::Vector{LineData}, 
        cch0, T::Float64, N::Float64, p::Float64)

    nb_species = length(linedata)

    lined    = Vector{Vector{SVector{9, Float64}}}(undef, nb_species)
    ΔλL_mean = Vector{Float64}(undef, nb_species)
    ΔλD_mean = Vector{Float64}(undef, nb_species)

    ispecies = 2
    for ispecies in eachindex(linedata)
        ldspecies, Qinpspecies, cspecies, mspecies = linedata[ispecies], moleculardata[ispecies].Qinp, cch0[ispecies], moleculardata[ispecies].iso_m;

        nλl = length(linedata[ispecies].λ210)
        ldl = Vector{SVector{9, Float64}}(undef, length(linedata[ispecies].λ210))
        ΔλLs = Vector{Float64}(undef, nλl)
        ΔλDs = Vector{Float64}(undef, nλl)

        iλ = 1
        #Threads.@threads 
        for iλ in 1:nλl
            mid, lid, λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ = 
                compute_line_emission_and_absorption_iλ(ldspecies, Qinpspecies, cspecies, mspecies, T, N, p, par.T_ref, iλ)

            ldl[iλ] = SVector{9,Float64}(λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ, mspecies[lid])
            ΔλLs[iλ] = ΔλL
            ΔλDs[iλ] = ΔλG
        end
        lined[ispecies] = ldl
        ΔλL_mean[ispecies] = Statistics.mean(ΔλLs)
        ΔλD_mean[ispecies] = Statistics.mean(ΔλDs)
    end

    lined, ΔλL_mean, ΔλD_mean
end

"""
    sum over all lines using their line shape

    T - temperature
    N - density
 """
function sum_over_lines(par, lined, T, λb)
    nbl = length(ldl)

    λ1 = λb[1]
    λend = λb[end]
    Δλ = λend - λ1
    dλ = λb[2] - λ1
    nλ = length(λb)

    κbt = zeros(Float64, nλ, Threads.nthreads())
    ϵbt = zeros(Float64, nλ, Threads.nthreads())
    fbt = zeros(Float64, nλ, Threads.nthreads())

    # λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ, mass, Float64(mid)

    iline = 1
    Threads.@threads for iline in eachindex(ldl)
        λ21  = lined[iline][1]      
        γ    = lined[iline][2]    
        ΔλL  = lined[iline][3]      
        ΔλG  = lined[iline][4]      
        N1   = lined[iline][5]     
        N2   = lined[iline][6]     
        ϵ    = lined[iline][7]    
        κ    = lined[iline][8]    
        mass = lined[iline][9]       

        if λ21 >= λ1 && λ21 <= λend
            iλ = floor(Int64, (λ21 - λ1) / Δλ * Float64(nλ-1)) + 1

            gauss   = GaussProfile(mass, T)
            lorentz = LorentzProfile(ΔλL)

            δiλ = max(2, floor(Int64, (ΔλL + ΔλG) * par.Δλ_factor / dλ))
            iλm = max(1, iλ - δiλ)
            iλp = min(iλ + δiλ + 1, nλ)

            sumft = zeros(Float64, Threads.threadid())
            for iλ in iλm:iλp
                fbt[iλ, Threads.threadid()] = voigt(gauss, lorentz, λb[iλ], λ21)
                sumft[Threads.threadid()] += fbt[iλ, Threads.threadid()]
            end
            int_fb = sum(sumft) * dλ
            cf = 1.0
            if int_fb > 0.8 && int_fb <= 1.0
                cf = 1.0/int_fb
            end

            for iλ in iλm:iλp
                κbt[iλ, Threads.threadid()] = κbt[iλ, Threads.threadid()] + κ * fbt[iλ, Threads.threadid()] * cf
                ϵbt[iλ, Threads.threadid()] = ϵbt[iλ, Threads.threadid()] + ϵ * fbt[iλ, Threads.threadid()] * cf
            end
        end
    end

    κb = zeros(Float64, nλ)
    ϵb = zeros(Float64, nλ)
    for tid in 1:Threads.nthreads()
        for iλ in 1:nλb
            κb[iλ] += κbt[iλ, tid]
            ϵb[iλ] += ϵbt[iλ, tid]
        end
    end

    κb, ϵb
end

function integrate_intensity_over_Δs(Iλb::Vector{Float64}, κb::Vector{Vector{Float64}}, ϵb::Vector{Vector{Float64}}, 
    Δs::Float64, with_emission::Bool, κΔs_limit::Float64)
    
    with_emission = par.with_emission
    iλ = 1
    for iλ in eachindex(Iλb)

        κbλ = 0.0
        ϵbλ = 0.0
        for ispecies in eachindex(κb)
            κbλ += κb[ispecies][iλ]
            ϵbλ += ϵb[ispecies][iλ]
        end

        exp_κ = exp(-κbλ * Δs)
        eps = 0.0
        if with_emission
            eps = 0.0
            if abs(κbλ) * Δs < κΔs_limit
                eps = ϵbλ * Δs
            else
                eps = ϵbλ / κbλ * (1.0 - exp_κ)
            end
            Iλb[iλ] = Iλb[iλ] * exp_κ + eps
        else
            Iλb[iλ] = Iλb[iλ] * exp_κ
        end

        if isnan(Iλb[iλ])
            @infoe @sprintf("%d  %e  %e  %e  %e", iλ, Iλb[iλ], κbsλ, ϵbsλ, exp_κ)
        end
    end
end

function add_background()
#        #@time begin
#        # add background ?
#        if spec.par.background > 1.0e-10
#            iw = floor(Int64, ((ΔλL_mean + ΔλD_mean) / dλ * spec.par.Δλ_factor))
#            # compute movning average
#            ma_κ = moving_average5(spec.κ_c, iw*2)
#            ma_ϵ = moving_average5(spec.ϵ_c, iw*2)
#            # add background
#            for iλ in 1:nb_λ
#                spec.κ_c[iλ] += ma_κ[iλ] * spec.par.background
#                spec.ϵ_c[iλ] += ma_ϵ[iλ] * spec.par.background
#            end
#            # ning average once
#            if iN == 1 && iθ == 1
#                save_intensity_as_hdf5(joinpath(spec.par.out_dir, "moving_average_kappa"), spec.λ, ma_κ)
#            end
#        end
end

function initial_intensity(par, λ)
    # initial intensity
    @infoe @sprintf("initial_intensity = %s", par.initial_intensity)
    Iλ = if par.initial_intensity == :planck
        Iλ = planck_λ(par.surface_T, λ)
        Iλ .* (1.0 - par.albedo)
    else
        zeros(Float64, length(λ))
    end
    Iλ
end

# integrate along the path
function integrate_along_path(par, atm, moleculardata, linedata, ch0, iθ, θ, λb)
    κΔs_limit = 0.01

    # number of lines wavelength intervall
    nλl = length(linedata[1].λ210)
    dλl = linedata[1].λ210[2] - linedata[1].λ210[1]

    # number of wavelengths and spacing
    nλb = length(λb)
    dλb = λb[2] - λb[1]

    # initial integrated intensity
    Iλb = initial_intensity(par, λb)
    int_I0 = sum(Iλb) * dλb
    #@infoe @sprintf("p = %e  T = %e  N = %e, int_I = %e, initial_intensity = %s", p, T, N, int_I0, par.initial_intensity)

    Tmin = par.surface_T
    Nmin = 1.0e30

    result_data = Results(19)

    ih = 1
    h = atm.h[ih]
    nh = length(atm.h)
    @time begin
    for ih in eachindex(atm.h)
        # >> 1
        cch0 = get_concentrations(atm, atm.h[ih], ch0)
        # << 1

        t1 = CPUtime_us()
        # pressure, temperature and density at height = z
        p = atm.p[ih]
        T = atm.T[ih]
        N = atm.N[ih]

        if par.T_of_h == false
            T = par.surface_T
            if par.N_of_h == false
                N = p / (c_kB * T)
            end
        end

        if T < Tmin
            Tmin = T
        end
        if N < Nmin
            Nmin = N
        end

        t2 = CPUtime_us()
    
        # compute the line coefficients
        # >> 2
        lined, ΔλL_mean, ΔλD_mean = compute_lines_emission_and_absorption(moleculardata, linedata, cch0, T, N, p);
        # << 2

        t3 = CPUtime_us()

        # >> 3
        κb = Vector{Vector{Float64}}(undef, length(lined))
        ϵb = Vector{Vector{Float64}}(undef, length(lined))
        # 3
        for ispecies in eachindex(lined)
            κb[ispecies], ϵb[ispecies] = sum_over_lines(par, lined[ispecies], T, λb)
        end
        # << 3

        t4 = CPUtime_us()
        length(κb)
        ll = @. ifelse(κb > 0.0, 1, 0)
        sum(ll)

        # 4
        add_background()
        t5 = CPUtime_us()

        # step size Δs = z/cos(θ)
        Δs = 1.0
        if ih < nh
            Δs = (atm.h[ih+1] - atm.h[ih]) / cos(θ)
        else
            Δs = (atm.h[ih] - atm.h[ih-1]) / cos(θ)
        end

        # 5
        integrate_intensity_over_Δs(Iλb, κb, ϵb, Δs, par.with_emission, κΔs_limit)
        #@infoe I_λ[div(nb_λ,2)], sum(I_λ)
        t6 = CPUtime_us()

        # save intensity and spectrum
        if par.z_iout[istep] == 1
            znext = par.h[istep] + Δz

            Iname = @sprintf("intensity_%03d_%03d_%d_%d_%5.3e", outid, istep, iN, iθ, par.z[istep])
            fname = joinpath(par.outdir, "intensity", Iname)
            save_intensity_as_hdf5(fname, λb, Iλb)

            cname = @sprintf("spectrum_%03d_%03d_%d_%d_%5.3e", outid, istep, iN, iθ, par.z[istep])
            fname = joinpath(par.outdir, "spectrum", cname)
            save_spectrum_as_hdf5(fname, par)
        end
        t7 = CPUtime_us()

        # add results
        if istep == 1
            int_I, int_ϵ, int_κ, int_Iκ = add_results(result_data, dλ, 0.0, NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, I_λ, vϵ, vκ, vIκ)
        end
        int_I, int_ϵ, int_κ, int_Iκ = add_results(result_data, dλ, par.h[istep], NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, I_λ, vϵ, vκ, vIκ)

        ## write result_data to log file
        out = @sprintf("iz = %3d, z = %12.5e,  NCO2 = %12.5e, θ = %12.5e, T = %12.5e, N = %12.5e, I = %12.5e, ϵ = %12.5e, κ = %12.5e, Iκ = %12.5e, ΔλL = %12.5e, ΔλD = %12.5e",
                                istep, spec.z[istep], NCO2*1.0e6, θ*180.0/π, T, N, int_I, int_ϵ, int_κ, int_Iκ, ΔλL_mean, ΔλD_mean)
        write(spec.logfile, out * "\n")
        flush(spec.logfile)
        @infoe out
        t8 = CPUtime_us()
        
        ts = [t1, t2, t3, t4, t5, t6, t7, t8]
        str = []
        for i in eachindex(ts[1:end-1])
            push!(str, @sprintf("%d:%8.2e", i, ts[i+1] - ts[i]))
        end
        @infoe join(str, ", ")
    end  # lop over z istep
    end # @time

    result_data
end

# Compute absorption
function integrate(par::RunParameter, paths::OutPaths, atm::Atmosphere, moleculardata::Vector{MolecularData}, linedata::Vector{LineData})
    mkpath(paths.outdir)
    mkpath(paths.intensity_dir)
    
    # compute Planck intensity
    λ1 = 1.0e-6
    λ2 = 1.0e-4
    nλ = 1000
    λP = collect(range(λ1, λ2, nλ))
    IP = planck_λ(par.surface_T, λP)
    save_planck_as_hdf5(joinpath(paths.intensity_dir, paths.planck_single), par.surface_T, λP, IP)

    nλb = floor(Int64, (par.λmax - par.λmin) / par.Δλb)
    λb = collect(range(par.λmin, par.λmax, nλb))
    IPb = compute_planck(par.planck_Ts, λb)
    save_planck_as_hdf5(joinpath(paths.intensity_dir, paths.planck_multi),  par.planck_Ts, λb, IPb)

    #logfile header
#    ncol = 11
#    out = @sprintf("nb_NCO2 = %d,  nb_θ = %d, nb_z = %d, nb_columns = %d", size(spec.NCO2,1), size(θ,1), size(spec.z,1), ncol)
#    write(spec.logfile, string(out, "\n"))
#    @infoe out

    # loop over CO2 concentrations
    outid = 0
    ic    = 1
    iθ, θ = 1, 0.0
    for ic in eachindex(par.c_ppm)
        ch0 = par.c_ppm[ic]
        # loop over angles
        for (iθ, θ) in enumerate(par.θ)
            outid = outid + 1

            # >> 
            fname = joinpath(spec.par.out_dir, @sprintf("result_%03d_%d_%d.hdf5", outid, ic, iθ))
            @infoe @sprintf("Results file : %s", fname)

            # intermediate log file header line
            out = @sprintf("# iN = %d, iθ = %d, NCO2 = %12.5e, θ = %12.5e", iN, iθ, NCO2*1.0e6, θ/π*180.0)
            write(spec.logfile, string(out, "\n"))
            @infoe out
            # <<

            # integrate along path
            @time result_data = integrate_along_path(par, atm, moleculardata, linedata, ch0, iθ, θ, λb)

            # save result_data
            write_result_data(result_data, fname)
            @infoe @sprintf("Results saved to %s", fname)
        end
    end
end
