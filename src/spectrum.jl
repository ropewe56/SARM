import Statistics
using StaticArrays
using CPUTime
using Interpolations

using PhysConst
using SimpleLog

"""
No5
"""
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
    lid   = ldspecies.lid[iλ]

    if lid > 11
        @warne mid, lid, λ210
    end
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
    compute_lines_emission_and_absorption(moleculardata, linedata, Nmolecules, T, N, p)

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

No4
"""
function compute_lines_emission_and_absorption(moleculardata::Vector{MolecularData}, linedata::Vector{LineData}, 
        cch0, T::Float64, N::Float64, p::Float64)

    nb_species = length(linedata)

    linedata_pTN = Vector{Vector{SVector{9, Float64}}}(undef, nb_species)
    ΔλL_mean = Vector{Float64}(undef, nb_species)
    ΔλD_mean = Vector{Float64}(undef, nb_species)

    ispecies = 2
    for ispecies in eachindex(linedata)
        ldspecies, Qinpspecies, cspecies, mspecies = linedata[ispecies], moleculardata[ispecies].Qinp, cch0[ispecies], moleculardata[ispecies].iso_m;

        nλl    = length(linedata[ispecies].λ210)
        ld_pTN = Vector{SVector{9, Float64}}(undef, length(linedata[ispecies].λ210))
        ΔλLs   = Vector{Float64}(undef, nλl)
        ΔλDs   = Vector{Float64}(undef, nλl)

        iλ = 1
        #Threads.@threads 
        for iλ in 1:nλl
            # No5
            mid, lid, λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ = 
                compute_line_emission_and_absorption_iλ(ldspecies, Qinpspecies, cspecies, mspecies, T, N, p, par.T_ref, iλ)

            ld_pTN[iλ] = SVector{9,Float64}(λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ, mspecies[lid])
            ΔλLs[iλ] = ΔλL
            ΔλDs[iλ] = ΔλG
        end
        linedata_pTN[ispecies] = ld_pTN
        ΔλL_mean[ispecies] = Statistics.mean(ΔλLs)
        ΔλD_mean[ispecies] = Statistics.mean(ΔλDs)
    end

    linedata_pTN, ΔλL_mean, ΔλD_mean
end

"""
    sum over all lines using their line shape

    T - temperature
    N - density
"""
function sum_over_lines(par, lines, T, λb)
    nbl = length(lines)

    λ1 = λb[1]
    λend = λb[end]
    Δλ  = λend - λ1
    dλ  = λb[2] - λ1
    nλb = length(λb)

    κbt = alloc2(par.prealloc, :κbt, nλb, Threads.nthreads(), true)
    ϵbt = alloc2(par.prealloc, :ϵbt, nλb, Threads.nthreads(), true)
    fbt = alloc2(par.prealloc, :fbt, nλb, Threads.nthreads(), true)

    # λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ, mass, Float64(mid)

    Threads.@threads for il in eachindex(lines)
        λ21  = lines[il][1]      
        γ    = lines[il][2]    
        ΔλL  = lines[il][3]      
        ΔλG  = lines[il][4]      
        N1   = lines[il][5]     
        N2   = lines[il][6]     
        ϵ    = lines[il][7]    
        κ    = lines[il][8]    
        mass = lines[il][9]       

        if λ21 >= λ1 && λ21 <= λend
            iλ = floor(Int64, (λ21 - λ1) / Δλ * Float64(nλb-1)) + 1

            gauss   = GaussProfile(mass, T)
            lorentz = LorentzProfile(ΔλL)

            δiλ = max(2, floor(Int64, (ΔλL + ΔλG) * par.Δλ_factor / dλ))
            iλm = max(1, iλ - δiλ)
            iλp = min(iλ + δiλ + 1, nλb)

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

    κb = alloc1(par.prealloc, :κb, nλb, true)
    ϵb = alloc1(par.prealloc, :ϵb, nλb, true)
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

"""
No3
"""
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

"""
    integrate_along_path(par, atm, moleculardata, linedata, ch0, ic, iθ, θ, λb)
No2
"""
function integrate_along_path(par, paths, atm, moleculardata, linedata, ch0, result_id, ic, iθ, θ, λb)

    # number of lines wavelength intervall
    nλl = length(linedata[1].λ210)
    dλl = linedata[1].λ210[2] - linedata[1].λ210[1]

    # number of wavelengths and spacing
    nλb = length(λb)

    # initial integrated intensity
    # No3
    Iλb = initial_intensity(par, λb)
    int_I0 = sum(Iλb) * par.Δλb

    Tmin = par.surface_T
    Nmin = 1.0e30

    results = Results(19)

    nh = length(atm.h)
    ih = 1
    
    logfio = open(paths.logfile, "w")

    @time begin

    for ih in eachindex(atm.h)
        t1 = CPUtime_us()

        # >> 1
        cch0 = get_concentrations(moleculardata, ch0)
        # << 1
        
        # >> 1  pressure, temperature and density at height = z
        p = atm.p[ih]
        T = atm.T[ih]
        N = atm.N[ih]

        if par.T_of_h == false 
            T = par.surface_T
            if par.N_of_h == false
                N = p / (c_kB * T)
            end
        end
        if T < Tmin Tmin = T end
        if N < Nmin Nmin = N end
        # << 1

        t2 = CPUtime_us()
        # >> 2
        # No4
        # linedata_pTN: Vector{Vector{SVector{9, Float64}}}(undef, nb_species)
        # SVector: λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ, mass
        linedata_pTN, ΔλL_mean, ΔλD_mean = compute_lines_emission_and_absorption(moleculardata, linedata, cch0, T, N, p);
        # << 2
        t3 = CPUtime_us()

        # >> No6
        κb = Vector{Vector{Float64}}(undef, length(linedata_pTN))
        ϵb = Vector{Vector{Float64}}(undef, length(linedata_pTN))
        for ispec in eachindex(linedata_pTN)
            lines = linedata_pTN[ispec]
            κb[ispec], ϵb[ispec] = sum_over_lines(par, lines, T, λb)
        end
        # << 3
        t4 = CPUtime_us()
        
        # 4
        add_background()
        t5 = CPUtime_us()

        # step size Δs = z/cos(θ)
        Δs = 1.0
        if ih < nh
            Δs = (atm.h[ih+1] - atm.h[ih]) / cos(θ)
        else
            par.κΔs_limit
            Δs = (atm.h[ih] - atm.h[ih-1]) / cos(θ)
        end

        # No7
        integrate_intensity_over_Δs(Iλb, κb, ϵb, Δs, par.with_emission, par.κΔs_limit) # 3
        t6 = CPUtime_us()

        if atm.h_iout[ih] == 1
            md = moleculardata[1]
            write_to_hdf5(paths, atm, result_id, ih, ic, iθ, md, λb, Iλb, κb, ϵb)
        end

        # add results
        # No8
        if ih == 1
            int_I, int_ϵ, int_κ, int_Iκ = add_results!(results, par, 0.0, cch0, θ, T, N, ΔλL_mean, ΔλD_mean, Iλb, ϵb, κb)
        end
        int_I, int_ϵ, int_κ, int_Iκ = add_results!(results, par, atm.h[ih], cch0, θ, T, N, ΔλL_mean, ΔλD_mean, Iλb, ϵb, κb)

        ## write results to log file
        local out
        for i in eachindex(cch0)
            out = @sprintf("ih = %3d, h = %12.5e,  c = %12.5e, θ = %12.5e, T = %12.5e, N = %12.5e, I = %12.5e, ϵ = %12.5e, κ = %12.5e, Iκ = %12.5e, ΔλL = %12.5e, ΔλD = %12.5e",
                                ih, atm.h[ih], cch0[i], θ*180.0/π, T, N, int_I, int_ϵ[i], int_κ[i], int_Iκ[i], ΔλL_mean[i], ΔλD_mean[i])
            write(logfio, out * "\n")
        end
        flush(logfio)
        @infoe out
        
        t7 = CPUtime_us()
        
        ts = [t1, t2, t3, t4, t5, t6, t7]
        str = ["cputime"]
        for i in eachindex(ts[1:end-1])
            push!(str, @sprintf("%d:%8.2e", i, ts[i+1] - ts[i]))
        end
        @infoe join(str, ", ")
    end  # lop over z ih
    
    end # @time

    results
end

"""
    function integrate(par::RunParameter, paths::OutPaths, atm::Atmosphere, moleculardata::Vector{MolecularData}, linedata::Vector{LineData})
No1
"""
function integrate(par::RunParameter, paths::OutPaths, atm::Atmosphere, moleculardata::Vector{MolecularData}, linedata::Vector{LineData})
    
    # compute Planck intensity
    λ1 = 1.0e-6
    λ2 = 1.0e-4
    nλ = 1000
    λP = collect(range(λ1, λ2, nλ))
    IP = planck_λ(par.surface_T, λP)
    save_planck_as_hdf5(joinpath(paths.intensity, paths.planck_single), par.surface_T, λP, IP)

    nλb = floor(Int64, (par.λmax - par.λmin) / par.Δλb)
    λb = collect(range(par.λmin, par.λmax, nλb))
    IPb = compute_planck(par.planck_Ts, λb)
    save_planck_as_hdf5(joinpath(paths.intensity, paths.planck_multi),  par.planck_Ts, λb, IPb)

    #logfile header
#    ncol = 11
#    out = @sprintf("nb_NCO2 = %d,  nb_θ = %d, nb_z = %d, nb_columns = %d", size(spec.NCO2,1), size(θ,1), size(spec.z,1), ncol)
#    write(spec.logfile, string(out, "\n"))
#    @infoe out

    # loop over CO2 concentrations
    result_id = 0
    ic    = 1
    iθ, θ = 1, 0.0
    for ic in 1:size(par.c_ppm,2)
        ch0 = par.c_ppm[:,ic]
        # loop over angles
        for (iθ, θ) in enumerate(par.θ)
            result_id = result_id + 1

            # >> 
            fname = joinpath(paths.results, @sprintf("result_%03d_%d_%d.hdf5", result_id, ic, iθ))
            @infoe @sprintf("Results file : %s", fname)

            # intermediate log file header line
            out = @sprintf("# ic = %d, iθ = %d, ch0 = %s, θ = %12.5e", ic, iθ, ch0, θ)
            write(paths.logfile, string(out, "\n"))
            @infoe out
            # <<

            # integrate along path
            # No2
            @time result_data = integrate_along_path(par, paths, atm, moleculardata, linedata, ch0, result_id, ic, iθ, θ, λb)

            # save result_data
            write_result_data(result_data, fname)
            @infoe @sprintf("Results saved to %s", fname)
        end
    end
end
