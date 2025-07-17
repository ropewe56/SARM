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
    sum over all lines using their line shape

    T - temperature
    N - density
    ML = ML[:CO2]
"""
function sum_over_lines(par, MLspec, T, λb)
    λ1   = λb[1]
    λend = λb[end]
    Δλ   = λend - λ1
    dλ   = λb[2] - λ1
    nλb  = length(λb)
    par.Δλ_factor = 10

    nbthreads = Threads.nthreads()
    κbt = alloc2(par.prealloc, :κbt, nλb, nbthreads, true)
    ϵbt = alloc2(par.prealloc, :ϵbt, nλb, nbthreads, true)
    κb = alloc1(par.prealloc, :κb, nλb, true)
    ϵb = alloc1(par.prealloc, :ϵb, nλb, true)

    λ21  = MLspec[3, :]
    index = @. ifelse(λ21 >= λ1 && λ21 <= λend, true, false)
    ML = MLspec[:,index]
    n1, nλl = size(ML)

    Threads.@threads for il in 1:nλl
        tid = Threads.threadid()

        λ21  = ML[ 3, il]
        ΔλL  = ML[ 5, il]
        ΔλG  = ML[ 6, il]
        mass = ML[ 9, il]
        ϵ    = ML[10, il]
        κ    = ML[11, il]

        iλb = floor(Int64, (λ21 - λ1) / Δλ * Float64(nλb-1)) + 1

        δiλ = max(2, floor(Int64, (ΔλL + ΔλG) * par.Δλ_factor / dλ))
        iλm = max(  1, iλb - δiλ)
        iλp = min(nλb, iλb + δiλ + 1)

        #fG = fgauss(λb[iλm:iλp], λb[iλb], ΔλG)
        #fL = florentz(λb[iλm:iλp], λb[iλb], ΔλL)
        
        #fb = florentz(λb[iλm:iλp], λb[iλb], ΔλL+ΔλG)
        fb = voigt(ΔλG, ΔλL, fg, fl, λ, λ0)

        κbt[iλm:iλp, tid] += @. κ * fb
        ϵbt[iλm:iλp, tid] += @. ϵ * fb
    end

    κb[:] = sum(κbt,dims=2)
    ϵb[:] = sum(ϵbt,dims=2)

    plt.plot(κb[:])
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
    cputimes = []
    nh = length(atm.h)

    ih = 1
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

        cihic = Dict{Symbol, Float64}()
        ML    = Dict{Symbol,Matrix{Float64}}()
        for (spec, cc) in par.c_ppm
            md   = moleculardata[spec]
            miso = md.iso_m       # Vector
            Qiso = md.Qisoh[:,ih] # Vector 
            Qref = md.Qref # Vector 
            cihic[spec] = md.cnh[ih] * cc[ic] * PPM

            nλl = length(linedata[spec].λ210)
            ML[spec] = Matrix{Float64}(undef, 12, nλl)
            compute_lines_emission_and_absorption!(ML[spec], par, linedata[spec], Qref, Qiso, miso, cihic[spec], T, N, p);
        end

        κbs      = Dict{Symbol, Vector{Float64}}()
        ϵbs      = Dict{Symbol, Vector{Float64}}()
        ΔλL_mean = Dict{Symbol, Float64}()
        ΔλD_mean = Dict{Symbol, Float64}()
        spec = :CO2
        cc = par.c_ppm[spec]
        MLspec = ML[spec]
        for (spec, cc) in par.c_ppm
            κbs[spec], ϵbs[spec] = sum_over_lines(par, ML[spec], T, λb)
            ΔλLs = ML[spec][3,:]
            ΔλDs = ML[spec][4,:]
            ΔλL_mean[spec] = Statistics.mean(ΔλLs)
            ΔλD_mean[spec] = Statistics.mean(ΔλDs)
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
            write_results_to_hdf5(par.paths, atm, ic, iθ, ih, ML, λb, Iλb, κb, ϵb, κbs, ϵbs)
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
function integrate(par::RunParameter, rdb, atm::Atmosphere, molecular_data::Dict{Symbol,MolecularData},  line_data::Dict{Symbol,LineData})
    nλb = floor(Int64, (par.λmax - par.λmin) / par.Δλb)
    λb  = collect(range(par.λmin, par.λmax, nλb))
    create_planck_spectrum(par, λb)

    cch0 = [par.c_ppm[k][1] for k in keys(par.c_ppm)]
    nbc = maximum([length(par.c_ppm[k]) for k in keys(par.c_ppm)])

    # loop over CO2 concentrations
    ic     = 1
    iθ, θ  = 1, 0.0
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

