import Statistics
using StaticArrays
using CPUTime
using Interpolations
using DataFrames

using PhysConst
using SimpleLog

function integrate_results(par::RunParameter, h, θ, T, N, ΔλL_mean, ΔλD_mean, Iλb, ϵb, κb)
    nλb = length(Iλb)

    nspec = length(keys(κb))

    int_ϵ  = Vector{Dict{Symbol,Float64}}(undef, 3)
    int_κ  = Vector{Dict{Symbol,Float64}}(undef, 3)
    int_Iκ = Vector{Dict{Symbol,Float64}}(undef, 3)
    for i in 1:3
        int_ϵ[i]  = Dict{Symbol,Float64}()
        int_κ[i]  = Dict{Symbol,Float64}()
        int_Iκ[i] = Dict{Symbol,Float64}()
    end

    # integrate over all wavelengths
    int_I = [sum(Iλb) * par.Δλb]

    Iκ = Dict{Symbol, Vector{Float64}}()
    for (i, spec) in enumerate(keys(κb))
        Iκ[spec] = Iλb .* κb[spec]
    end

    for (i, spec) in enumerate(keys(κb))
        int_ϵ[1][spec]  = sum(ϵb[spec]) * par.Δλb
        int_Iκ[1][spec] = sum(Iκ[spec]) * par.Δλb
        int_κ[1][spec] = Statistics.mean(κb[spec])
    end

    # integrate over nλb/6,...,nλb-nλb/6  wavelengths
    n1 = floor(Int64, nλb/6)
    n2 = nλb-n1
    push!(int_I, sum(Iλb[n1:n2]) * par.Δλb)
    for (i, spec) in enumerate(keys(κb))
        int_ϵ[2][spec]  = sum(ϵb[spec][n1:n2]) * par.Δλb
        int_Iκ[2][spec] = sum(Iκ[spec][n1:n2]) * par.Δλb
        int_κ[2][spec]  = Statistics.mean(κb[spec][n1:n2])
    end

    # integrate over nλb/4,...,nλb-nλb/4  wavelengths
    n1 = floor(Int64, nλb/4)
    n2 = nλb-n1
    push!(int_I, sum(Iλb[n1:n2]) * par.Δλb)
    for (i, spec) in enumerate(keys(κb))
        int_ϵ[3][spec]  = sum(ϵb[spec][n1:n2]) * par.Δλb
        int_Iκ[3][spec] = sum(Iκ[spec][n1:n2]) * par.Δλb
        int_κ[3][spec]  = Statistics.mean(κb[spec][n1:n2])
    end

    int_I, int_ϵ, int_κ, int_Iκ
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
No5
"""
function compute_line_emission_and_absorption_iλ(ld::LineData, Qiso, miso, c, T, N, p, T_ref, iλ)
    dΩ = 1.0
    β  = 1.0/(c_kB * T)

    #                  1      2   3     4  5     6    7    8    9    10   11   12    13    14    15     16
    #lines[i] = SA_F64[λ_ul0, E_l, E_u, S, A_ul, γ_a, γ_s, n_a, δ_a, g_u, g_l, B_ul, B_lu, ΔλL0, iso_m, iso_c]

    λ210  = ld.λ210[iλ]
    E1    = ld.E1[iλ]
    E2    = ld.E2[iλ]
    A21   = ld.A21[iλ]
    γair  = ld.γair[iλ]
    γself = ld.γself[iλ]
    nair  = ld.nair[iλ]
    δair  = ld.δair[iλ]
    g2    = ld.g2[iλ]
    g1    = ld.g1[iλ]
    B12   = ld.B12[iλ]
    B21   = ld.B12[iλ]
    iso   = ld.lid[iλ]

    if iso > 11
        @warne mid, lid, λ210
    end
    λ21 = λ210 / (1.0 + λ210 * δair * p)

    # γ = (Tref/T)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})
    dT = (T_ref/T)^nair
    γ = dT * (γair * p * (1.0 - c) + γself * p * c)

    ΔλL = λ21^2 * γ
    ΔλG = sqrt(2.0 * c_kB * T / miso[iso]) / c_c * λ21

    N1  = g1 * exp(- E1 * β) / Qiso[iso] * c * N
    N2  = g2 * exp(- E2 * β) / Qiso[iso] * c * N

    # ϵ_λ = h * c / λ_0 / (4 * π) * N_u * A_ul * f_λ
    ϵ = c_h * c_c / λ21 *  N2 * A21 * dΩ / (4.0 * π)

    # κ_λ = h / λ_0 * N_l * B_lu * (1 - N_u/N_l * g_l/g_u) * λ_0**2 / c * f_λ
    κ = c_h * λ21 / c_c * (N1 * B12 - N2 * B21)

    iso, λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ
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
function compute_lines_emission_and_absorption(par, ld::LineData, Qiso, miso, c, T, N, p)
    nλl = length(ld.λ210)

    ld_pTN = Vector{SVector{9, Float64}}(undef, nλl)
    ΔλLs = Vector{Float64}(undef, nλl)
    ΔλDs = Vector{Float64}(undef, nλl)

    iλ = 1
    Threads.@threads for iλ in 1:nλl
        iso, λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ = compute_line_emission_and_absorption_iλ(ld, Qiso, miso, c, T, N, p, par.T_ref, iλ)
        ld_pTN[iλ] = SVector{9,Float64}(λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ, miso[iso])
        ΔλLs[iλ] = ΔλL
        ΔλDs[iλ] = ΔλG
    end
    ld_pTN, ΔλLs, ΔλDs
end

"""
    sum over all lines using their line shape

    T - temperature
    N - density
"""
function sum_over_lines(par, lines, T, λb)
    λ1   = λb[1]
    λend = λb[end]
    Δλ   = λend - λ1
    dλ   = λb[2] - λ1
    nλb  = length(λb)

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

function integrate_intensity_over_Δs(Iλb::Vector{Float64}, κbs::Dict{Symbol,Vector{Float64}}, ϵbs::Dict{Symbol,Vector{Float64}}, 
                                     Δs::Float64, with_emission::Bool, κΔs_limit::Float64)

    nλb = length(Iλb)
    κb = zeros(Float64, nλb)
    ϵb = zeros(Float64, nλb)
    for (k, val) in κbs
        @. κb += val
    end
    for (k, val) in ϵbs
        @. ϵb += val
    end

    Threads.@threads for iλ in eachindex(Iλb)
        exp_κ = exp(-κb[iλ] * Δs)
        eps = 0.0
        if with_emission
            eps = if abs(κb[iλ]) * Δs < κΔs_limit
                ϵb[iλ] * Δs
            else
                ϵb[iλ] / κb[iλ] * (1.0 - exp_κ)
            end
            Iλb[iλ] = Iλb[iλ] * exp_κ + eps
        else
            Iλb[iλ] = Iλb[iλ] * exp_κ
        end

        if isnan(Iλb[iλ])
            @infoe @sprintf("%d  %e  %e  %e  %e", iλ, Iλb[iλ], κb[iλ], ϵb[iλ], exp_κ)
        end
    end
    κb, ϵb
end

"""
    integrate_along_path(par, atm, moleculardata, linedata, ch0, ic, iθ, θ, λb)
No2
"""
function integrate_along_path(par, rdb, atm, moleculardata, linedata, ic, iθ, θ, λb)
    nλb = length(λb)
    Iλb = initial_intensity(par, λb)
    int_I0 = sum(Iλb) * par.Δλb

    Tmin = par.surface_T
    Nmin = 1.0e30

    logfio = open(par.paths.logfile, "w")
    ih = 1
    h = atm.h[ih]
    cputimes = []
    nh = length(atm.h)
    for (ih,h) in enumerate(atm.h)
        tt = [time_ns()]
        # >> 1  pressure, temperature and density at height = z
        p = atm.p[ih]
        T = atm.T[ih]
        N = atm.N[ih]
        if par.T_of_h == false 
            Th = par.surface_T
            if par.N_of_h == false
                N = p / (c_kB * T)
            end
        end
        if T < Tmin Tmin = T end
        if N < Nmin Nmin = N end
        # << 1

        push!(tt, time_ns())
        # >> 2
        # No4
        # linedata_pTN: Vector{Vector{SVector{9, Float64}}}(undef, nb_species)
        # SVector: λ21, γ, ΔλL, ΔλG, N1, N2, ϵ, κ, mass
        κbs = Dict{Symbol, Vector{Float64}}()
        ϵbs = Dict{Symbol, Vector{Float64}}()
        ΔλL_mean = Dict{Symbol, Float64}()
        ΔλD_mean = Dict{Symbol, Float64}()
        cihic = Dict{Symbol, Float64}()
        for (spec, cc) in par.c_ppm
            md   = moleculardata[spec]
            miso = md.iso_m  # Vector
            Qiso = md.Qisoh[:,ih] # Vector 
            c    = md.cnh[ih] * cc[ic] * PPM
            ld   = linedata[spec]

            ld_pTN, ΔλLs, ΔλDs = compute_lines_emission_and_absorption(par, ld, Qiso, miso, c, T, N, p);
            κbs[spec], ϵbs[spec] = sum_over_lines(par, ld_pTN, T, λb)
            ΔλL_mean[spec] = Statistics.mean(ΔλLs)
            ΔλD_mean[spec] = Statistics.mean(ΔλDs)
            cihic[spec] = c
        end
        # << 2
        
        push!(tt, time_ns())
        
        ## 4
        #add_background()
        #t5 = time_ns()

        # step size Δs = z/cos(θ)
        Δs = 1.0
        if ih < nh
            Δs = (atm.h[ih+1] - atm.h[ih]) / cos(θ)
        else
            par.κΔs_limit
            Δs = (atm.h[ih] - atm.h[ih-1]) / cos(θ)
        end

        # No7
        κb, ϵb = integrate_intensity_over_Δs(Iλb, κbs, ϵbs, Δs, par.with_emission, par.κΔs_limit) # 3
        push!(tt, time_ns())

        hdf5_path = if atm.h_iout[ih] == 1
            write_results_to_hdf5(par.paths, atm, ic, iθ, ih, λb, Iλb, κb, ϵb, κbs, ϵbs)
        else
            missing
        end
        push!(tt, time_ns())

        # add results
        int_I, int_ϵ, int_κ, int_Iκ = integrate_results(par, atm.h[ih], θ, T, N, ΔλL_mean, ΔλD_mean, Iλb, ϵbs, κbs)        
        insert_into_rdb(rdb, ic, iθ, ih, atm.h[ih], θ, T, N, cihic, ΔλL_mean, ΔλD_mean, int_I, int_ϵ, int_κ, int_Iκ, hdf5_path)


        ## write results to log file
        for (i, spec) in enumerate(keys(cihic))
            out = @sprintf("%s, ih = %3d, h = %12.5e,  c = %12.5e, θ = %12.5e, T = %12.5e, N = %12.5e, I = %12.5e, ϵ = %12.5e, κ = %12.5e, Iκ = %12.5e, ΔλL = %12.5e, ΔλD = %12.5e",
                                spec, ih, atm.h[ih], cihic[spec], θ*180.0/π, T, N, int_I[1], int_ϵ[1][spec], int_κ[1][spec], int_Iκ[1][spec], ΔλL_mean[spec], ΔλD_mean[spec])
            @infoe out
            write(logfio, out * "\n")
        end
        flush(logfio)
        
        push!(tt, time_ns())

        dt = tt[2:end] - tt[1:end-1]
        push!(cputimes, [Float64(x).*1.0e-6 for x in dt])
    end  # lop over z ih
    reduce(hcat, cputimes)'
end

"""
    function integrate(par::RunParameter, atm::Atmosphere, moleculardata::Vector{MolecularData}, linedata::Vector{LineData})
No1
"""
function integrate(par::RunParameter, atm::Atmosphere, molecular_data::Dict{Symbol,MolecularData},  line_data::Dict{Symbol,LineData})
    nλb = floor(Int64, (par.λmax - par.λmin) / par.Δλb)
    λb  = collect(range(par.λmin, par.λmax, nλb))
    create_planck_spectrum(par, λb)

    # loop over CO2 concentrations
    ic     = 1
    iθ, θ  = 1, 0.0

    cch0 = [par.c_ppm[k][1] for k in keys(par.c_ppm)]
    nbc = maximum([length(par.c_ppm[k]) for k in keys(par.c_ppm)])
    rdb = ResultDB(par.paths.dbpath, par.c_ppm)
    if rdb === nothing
        return
    end

    for ic in 1:nbc
        # loop over angles
        for (iθ, θ) in enumerate(par.θ)
            # >> 
            out = @sprintf("# ic = %d, iθ = %d, ch0 = %s, θ = %12.5e", ic, iθ, cch0, θ)
            write(par.paths.logfile, string(out, "\n"))
            @infoe out
            # <<

            # integrate along path
            @time cputimes = integrate_along_path(par, rdb, atm, moleculardata, linedata, ic, iθ, θ, λb);            

            m1, m2 = size(cputimes)
            for im in 1:m2
                tim = cputimes[:,im]
                @printf("%d : %8.2e ms, %8.2e ms\n", im, sum(tim), Statistics.mean(tim))
            end 
            @printf("sum: %8.2e ms\n", sum(cputimes))
        end
    end
end

