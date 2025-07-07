import Statistics
using StaticArrays
using CPUTime

using PhysConst
using SimpleLog

mutable struct Spectrum
    line_data   :: LineData
    T_ref       :: Float64          # reference temperature
    p_ref       :: Float64          # reference pressure
    Q_CO2       :: Interpolator     # the Interplator for the
    p_vs_h      :: Interpolator     # pressure over height
    t_vs_h      :: Interpolator     # temperature over height
    NCO2        :: Vector{Float64}  # CO2 concentration values to iterate over
    θdeg        :: Vector{Float64}  # the angles to iterate over
    z           :: Vector{Float64}  # height values
    z_iout      :: Vector{Int64}    # where to save intensites
    λ           :: Vector{Float64}  # wavelengths
    κ_c         :: Vector{Float64}  # absorptions over wavelengths
    ϵ_c         :: Vector{Float64}  # emissivity  over wavelengths
    par         :: Parameter        # Parameter object
    logfile     :: IOStream         # logfile
end

function Spectrum(parameter_json_file)
    # create a Parameter object form the jso strinh
    par = Parameter(parameter_json_file)

    line_data = LineData(par.specdat_path, par.max_isotope_id)

    # height values
    z = read_npy(par.z_path)
    z_iout = read_npy(par.z_iout)

    # CO2 concentration values
    NCO2 = read_npy(par.NCO2_path)

    # angles
    θdeg = read_npy(par.theta_path)

    # CO2 partition functions
    Q_CO2 = Interpolator(par.T_Q_path)

    # pressure over height
    p_vs_h = Interpolator(par.h_p_path)

    # temperature over height
    t_vs_h = Interpolator(par.h_T_path)

    # wavelengths array
    nb_λ = floor(Int64, ((par.λmax - par.λmin) / par.Δλ))
    λ = collect(LinRange(par.λmin, par.λmax, nb_λ))

    # logfile
    @infoe joinpath(par.out_dir, par.logfile)
    logfile = open(joinpath(par.out_dir, par.logfile), "w")

    κ_c = zeros(Float64, nb_λ)
    ϵ_c = zeros(Float64, nb_λ)

    Spectrum(   line_data  ,
                par.T_ref  ,
                par.p_ref  ,
                Q_CO2      ,
                p_vs_h     ,
                t_vs_h     ,
                NCO2       ,
                θdeg       ,
                z          ,
                z_iout     ,
                λ          ,
                κ_c        ,
                ϵ_c        ,
                par        ,
                logfile    )
end

# Save convolved absorption and emission spectra to hdf5 file
function save_convolved_to_hdf5(spec::Spectrum, path)
    nb_λ = size(spec.λ, 1)
    a = Matrix{Float64}(undef, 3, nb_λ)
    for i in 1:nb_λ
        a[1, i] = spec.λ[i]
        a[2, i] = spec.κ_c[i]
        a[3, i] = spec.ϵ_c[i]
    end
    hdf5_path = string(path, ".hdf5")
    save_array_as_hdf5(hdf5_path, a, script_dir=false)
end

# Save convolved absorption and emission spectra to ny file
function save_spectrum_as_npy(path, spec::Spectrum)
    nb_λ = size(spec.λ, 1)
    a = Matrix{Float64}(undef, 3, nb_λ)
    for i in 1:nb_λ
        a[1,i] = spec.λ[i];
        a[2,i] = spec.κ_c[i];
        a[3,i] = spec.ϵ_c[i];
    end
    npy_path = string(path, ".hdf5")
    write_npy(npy_path, a)
end

function save_spectrum_as_hdf5(path, spec::Spectrum)
    nb_λ = size(spec.λ, 1)
    A = Matrix{Float64}(undef, 3, nb_λ)
    for i in 1:nb_λ
        A[1,i] = spec.λ[i];
        A[2,i] = spec.κ_c[i];
        A[3,i] = spec.ϵ_c[i];
    end
    hdf5_path = string(path, ".hdf5")
    save_array_as_hdf5(hdf5_path, A, script_dir=false)
end

function test_intepolator(spec::Spectrum)
    nb_zsteps = size(spec.z,1)
    for istep in 1:nb_zsteps
        p = y_at(spec.p_vs_h, spec.z[istep])
        T = y_at(spec.t_vs_h, spec.z[istep])
        N = p / (c_kB * T)
        @info @sprintf("z = %8.2e, p = %8.2e, T = %8.2e, N = %8.2e", spec.z[istep], p, T, N)
    end
end


function compute_emission_and_absorption_(line_data::LineData, Q_CO2, NCO2, T, N, p, T_ref, iline)
    ldc = line_data.linesc

    dΩ = 1.0
    β  = 1.0/(c_kB * T)

    #                  1      2   3     4  5     6    7    8    9    10   11   12    13    14    15     16
    #lines[i] = SA_F64[λ_ul0, E_l, E_u, S, A_ul, γ_a, γ_s, n_a, δ_a, g_u, g_l, B_ul, B_lu, ΔλL0, iso_m, iso_c]

    λ_ul0  = ldc[i_λ_ul0, iline]
    E_l    = ldc[i_E_l  , iline]
    E_u    = ldc[i_E_u  , iline]
    A_ul   = ldc[i_A_ul , iline]
    γ_a    = ldc[i_γ_a  , iline]
    γ_s    = ldc[i_γ_s  , iline]
    n_a    = ldc[i_n_a  , iline]
    δ_a    = ldc[i_δ_a  , iline]
    g_u    = ldc[i_g_u  , iline]
    g_l    = ldc[i_g_l  , iline]
    iso_m  = ldc[i_iso_m, iline]
    iso_c  = ldc[i_iso_c, iline]
    B_lu   = ldc[i_B_lu , iline]
    B_ul   = ldc[i_B_ul , iline]

    iso_id = line_data.iso_ids[iline]
    Q  = y_at(Q_CO2, T, ir=iso_id)

    Niso = N * NCO2 * iso_c

    λ_ul = λ_ul0 / (1.0 + λ_ul0 * δ_a * p)

    # γ = (Tref/T)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})
    dT = (T_ref/T)^n_a
    γ = dT * (γ_a * p * (1.0 - NCO2) + γ_s * p * NCO2)

    ΔλL = λ_ul^2 * γ
    ΔλG = sqrt(2.0 * c_kB * T / iso_m) / c_c * λ_ul

    N_l  = g_l * exp(- E_l * β) / Q * Niso
    N_u  = g_u * exp(- E_u * β) / Q * Niso

    # ϵ_λ = h * c / λ_0 / (4 * π) * N_u * A_ul * f_λ
    ϵ = c_h * c_c / λ_ul *  N_u * A_ul * dΩ / (4.0 * π)

    # κ_λ = h / λ_0 * N_l * B_lu * (1 - N_u/N_l * g_l/g_u) * λ_0**2 / c * f_λ
    κ = c_h * λ_ul / c_c * (N_l * B_lu - N_u * B_ul)

    ldv = line_data.linesv
    ldv[i_γ   ,iline] = γ
    ldv[i_ΔλL ,iline] = ΔλL
    ldv[i_ΔλG ,iline] = ΔλG
    ldv[i_N_l ,iline] = N_l
    ldv[i_N_u ,iline] = N_u
    ldv[i_ϵ   ,iline] = ϵ
    ldv[i_κ   ,iline] = κ
    ldv[i_λ_ul,iline] = λ_ul

    ΔλL, ΔλG
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
dT = (spec.par.T_ref/T)^n_a
γ = dT * (γ_a * p * (1.0 - NCO2) + γ_s * p * NCO2)

ΔλL = λ_ul^2 * γ
ΔλG = sqrt(2.0 * c_kB * T / iso_m) / c_c * λ_ul
ΔλL_mean[iline] = ΔλL
ΔλD_mean[iline] = ΔλG

$N_l  = \dfrac{g_l}{Q(T, iso)} \exp(- E_l  β)  N_{iso}$
$N_u  = \dfrac{g_u}{Q(T, iso)} \exp(- E_u  β)  N_{iso}$

$ϵ(λ) = \dfrac{h c}{λ_0} N_u A_{ul} * f(λ) \dfrac{dΩ}{4 π}$
$κ(λ) = \dfrac{h c}{λ_0} N_l B_{lu}  \left(1 - \dfrac{N_u}{N_l}  \dfrac{g_l}{g_u}\right)  \dfrac{λ_0^2}{c} f(λ)$

"""
function compute_emission_and_absorption(spec::Spectrum, T::Float64, N::Float64, p::Float64, NCO2::Float64)
    flag = false
    nb_lines = spec.line_data.nb_lines
    ΔλL_mean = Vector{Float64}(undef, nb_lines)
    ΔλD_mean = Vector{Float64}(undef, nb_lines)

    Threads.@threads for iline in 1:nb_lines
        ΔλL_mean[iline], ΔλD_mean[iline] = compute_emission_and_absorption_(spec.line_data, spec.Q_CO2, NCO2, T, N, p, spec.par.T_ref, iline)
    end

    # return average line widths
    Statistics.mean(ΔλL_mean), Statistics.mean(ΔλD_mean)
end

"""
    sum over all lines using their line shape

    T - temperature
    N - density
 """
function sum_over_lines(spec::Spectrum, T, N)
    nb_lines = spec.line_data.nb_lines
    nb_λ = size(spec.λ,1)

    λ_1 = spec.λ[1]
    λ_2 = spec.λ[nb_λ]
    Δλ  = spec.λ[nb_λ] - λ_1
    dλ  = spec.λ[2] - λ_1

    κ_c = zeros(Float64, nb_λ, Threads.nthreads())
    ϵ_c = zeros(Float64, nb_λ, Threads.nthreads())
    ff  = zeros(Float64, nb_λ, Threads.nthreads())

    ldc = spec.line_data.linesc
    ldv = spec.line_data.linesv

    # sum over lines
    Threads.@threads for iline in 1:nb_lines

        #iso_id = spec.line_data.iso_ids[iline]

        λ_ul = ldv[i_λ_ul, iline]

        if λ_ul >= λ_1 && λ_ul <= λ_2
            iλ0 = floor(Int64, (λ_ul - λ_1) / Δλ * Float64(nb_λ-1) )

            gauss = GaussProfile(ldc[i_iso_m, iline], T, λ_ul)
            lorentz = LorentzProfile(ldv[i_ΔλL, iline])

            diλ = max(2, floor(Int64, (ldv[i_ΔλL, iline] + ldv[i_ΔλG, iline]) * spec.par.Δλ_factor / dλ))
            iλm = max(1, iλ0 - diλ)
            iλp = min(iλ0 + diλ + 1, nb_λ-1)

            κ = ldv[i_κ, iline]
            ϵ = ldv[i_ϵ, iline]

            int_f = 0.0
            for iλ in iλm:iλp
                ff[iλ, Threads.threadid()] = voigt(gauss, lorentz, spec.λ[iλ], λ_ul)
                int_f += ff[iλ, Threads.threadid()]
            end
            int_f *= dλ
            cf = 1.0
            if int_f > 0.8 && int_f <= 1.0
                cf = 1.0/int_f
            end

            for iλ in iλm:iλp
                κ_c[iλ, Threads.threadid()] = κ_c[iλ, Threads.threadid()] + κ * ff[iλ, Threads.threadid()] * cf
                ϵ_c[iλ, Threads.threadid()] = ϵ_c[iλ, Threads.threadid()] + ϵ * ff[iλ, Threads.threadid()] * cf
            end
        end
    end

    spec.κ_c = zeros(Float64, nb_λ)
    spec.ϵ_c = zeros(Float64, nb_λ)
    for tid in 1:Threads.nthreads()
        for iλ in 1:nb_λ
            spec.κ_c[iλ] += κ_c[iλ, tid]
            spec.ϵ_c[iλ] += ϵ_c[iλ, tid]
        end
    end

end

function integrate_over_lines(nb_λ::Int64, I_λ::Vector{Float64}, κ_c::Vector{Float64}, ϵ_c::Vector{Float64}, Δs::Float64, with_emission_p::Bool, κds_limit::Float64)
    # integrate from s => s+Δs at all wavelengths
    vIκ = Vector{Float64}(undef, nb_λ)
    vκ  = Vector{Float64}(undef, nb_λ)
    vϵ  = Vector{Float64}(undef, nb_λ)
    for iλ in 1:nb_λ
        κ = κ_c[iλ]
        ϵ = ϵ_c[iλ]
        #tot_e = tot_e + ϵ * Δs
        exp_κ = exp(-κ * Δs)

        if with_emission_p
            eps = 0.0
            if abs(κ) * Δs < κds_limit
                eps = ϵ * Δs
            else
                eps = ϵ / κ * (1.0 - exp_κ)
            end
            I_λ[iλ] = I_λ[iλ] * exp_κ + eps
        else
            I_λ[iλ] = I_λ[iλ] * exp_κ
        end

        if isnan(I_λ[iλ])
            @infoe @sprintf("%d  %e  %e  %e  %e", iλ, I_λ[iλ], κ, ϵ, exp_κ)
        end

        vϵ[iλ]  = ϵ
        vκ[iλ]  = κ
        vIκ[iλ] = I_λ[iλ] * κ
    end # loop over λ
    vϵ, vκ, vIκ
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

function initial_intensity(spec::Spectrum)
    # initial intensity
    @infoe @sprintf("initial_intensity = %s", spec.par.initial_intensity)
    if spec.par.initial_intensity == "planck"
        I_λ = planck_λ(spec.par.T_surface, spec.λ)
        I_λ .* (1.0 - spec.par.albedo)
    else
        zeros(Float64, length(spec.λ))
    end
end

# integrate along the path
function integrate_along_path(spec::Spectrum, outid, iN, iθ, NCO2, θ)
    nb_zsteps = size(spec.z,1)
    Δz = 2.0e3
    κds_limit = 0.01

    # number of wavelengths and wavelength intervall
    nb_λ = length(spec.λ)
    dλ = spec.λ[2] - spec.λ[1]

    I_λ = initial_intensity(spec)

    # pessure, temperature, density at z = 0
    p = y_at(spec.p_vs_h, 0.0)
    T = y_at(spec.t_vs_h, 0.0)
    N = p / (c_kB * T)

    # initial integrated intensity
    int_I0 = sum(I_λ) * dλ
    @infoe @sprintf("p = %e  T = %e  N = %e, int_I = %e, initial_intensity = %s", p, T, N, int_I0, spec.par.initial_intensity)

    Tmin = spec.par.T_surface
    Nmin = 1.0e30

    result_data = Results(19)

    @time begin
    for istep in 1:nb_zsteps
        t1 = CPUtime_us()
        # pressure, temperature and density at height = z
        p = y_at(spec.p_vs_h, spec.z[istep])
        T = y_at(spec.t_vs_h, spec.z[istep])
        N = p / (c_kB * T)

        if spec.par.T_of_h == false
            T = spec.par.T_surface
            if spec.par.N_of_h == false
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
        (ΔλL_mean, ΔλD_mean) = compute_emission_and_absorption(spec, T, N, p, NCO2)
        t3 = CPUtime_us()

        sum_over_lines(spec, T, N)
        t4 = CPUtime_us()

        add_background()
        t5 = CPUtime_us()

        # step size Δs = z/cos(θ)
        Δs = 1.0
        if istep < nb_zsteps
            Δs = (spec.z[istep+1] - spec.z[istep]) / cos(θ)
        else
            Δs = (spec.z[istep] - spec.z[istep-1]) / cos(θ)
        end

        vϵ, vκ, vIκ = integrate_over_lines(nb_λ, I_λ, spec.κ_c, spec.ϵ_c, Δs, spec.par.with_emission, κds_limit)
        #@infoe I_λ[div(nb_λ,2)], sum(I_λ)
        t6 = CPUtime_us()

        # save intensity and spectrum
        if spec.z_iout[istep] == 1
            znext = spec.z[istep] + Δz

            Iname = @sprintf("intensity_%03d_%03d_%d_%d_%5.3e", outid, istep, iN, iθ, spec.z[istep])
            fname = joinpath(spec.par.out_dir, "intensity", Iname)
            save_intensity_as_hdf5(fname, spec.λ, I_λ)

            cname = @sprintf("spectrum_%03d_%03d_%d_%d_%5.3e", outid, istep, iN, iθ, spec.z[istep])
            fname = joinpath(spec.par.out_dir, "spectrum", cname)
            save_spectrum_as_hdf5(fname, spec)
        end
        t7 = CPUtime_us()

        # add results
        if istep == 1
            int_I, int_ϵ, int_κ, int_Iκ = add_results(result_data, dλ, 0.0, NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, I_λ, vϵ, vκ, vIκ)
        end
        int_I, int_ϵ, int_κ, int_Iκ = add_results(result_data, dλ, spec.z[istep], NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, I_λ, vϵ, vκ, vIκ)

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
function integrate(spec::Spectrum)
    # compute Planck intensity
    λ1 = 1.0e-6
    λ2 = 1.0e-4
    nλ = 1000
    T = spec.par.T_surface
    mkpath(spec.par.out_dir)
    compute_and_save_planck(joinpath(spec.par.out_dir, "planck_intensity"), λ1, λ2, nλ, T)
    compute_and_save_planck(joinpath(spec.par.out_dir, "planck_intensity_part"), spec.λ[1], spec.λ[end], size(spec.λ,1), spec.par.Planck_Ts)

    # vector of angles
    vθ = spec.θdeg * π/180.
    # vector of CO2 concentrations
    vNCO2 = copy(spec.NCO2)

    #logfile header
    ncol = 11
    out = @sprintf("nb_NCO2 = %d,  nb_θ = %d, nb_z = %d, nb_columns = %d", size(spec.NCO2,1), size(vθ,1), size(spec.z,1), ncol)
    write(spec.logfile, string(out, "\n"))
    @infoe out

    # loop over CO2 concentrations
    outid = 0
    for (iN, NCO2) in enumerate(vNCO2)
        # loop over angles
        for (iθ, θ) in enumerate(vθ)
            outid = outid + 1

            fname = joinpath(spec.par.out_dir, @sprintf("result_%03d_%d_%d.hdf5", outid, iN, iθ))
            @infoe @sprintf("Results file : %s", fname)

            # intermediate log file header line
            out = @sprintf("# iN = %d, iθ = %d, NCO2 = %12.5e, θ = %12.5e", iN, iθ, NCO2*1.0e6, θ/π*180.0)
            write(spec.logfile, string(out, "\n"))
            @infoe out

            # integrate along path
            @time result_data = integrate_along_path(spec, outid, iN, iθ, NCO2, θ)

            # save result_data
            write_result_data(result_data, fname)
            @infoe @sprintf("Results saved to %s", fname)
        end
    end
end
