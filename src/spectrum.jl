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


"""
    integrate_along_path(par, atm, moleculardata, linedata, ch0, ic, iθ, θ)
No2
"""
function integrate_along_path(par, rdb, atm, moleculardata, linedata, ic, iθ, θ)

    Iλb = initial_intensity(par)
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
        linedata_pTNc = Dict{Symbol,Matrix{Float64}}()
        for (spec, cc) in par.c_ppm
            md   = moleculardata[spec]
            miso = md.iso_m       # Vector
            Qiso = md.Qisoh[:,ih] # Vector 
            Qref = md.Qref # Vector 
            cihic[spec] = md.cnh[ih] * cc[ic] * PPM

            nλl = length(linedata[spec].λ210)
            linedata_pTNc[spec] = Matrix{Float64}(undef, 12, nλl)
            compute_lines_emission_and_absorption!(linedata_pTNc[spec], par, linedata[spec], Qref, Qiso, miso, cihic[spec], T, N, p);
        end

        κbs      = Dict{Symbol, Vector{Float64}}()
        ϵbs      = Dict{Symbol, Vector{Float64}}()
        ΔλL_mean = Dict{Symbol, Float64}()
        ΔλD_mean = Dict{Symbol, Float64}()
        spec = :CO2
        cc = par.c_ppm[spec]
        linedata_pTNc_spec = linedata_pTNc[spec]
        for (spec, cc) in par.c_ppm
            κbs[spec], ϵbs[spec] = sum_over_lines(par, linedata_pTNc[spec])
            ΔλLs = linedata_pTNc[spec][3,:]
            ΔλDs = linedata_pTNc[spec][4,:]
            ΔλL_mean[spec] = Statistics.mean(ΔλLs)
            ΔλD_mean[spec] = Statistics.mean(ΔλDs)
        end
        # << 2
        
        push!(tt, time_ns())
        
        ## 4
        #add_background()
        #t5 = time_ns()

        # step size Δs = z/cos(θ)
        Δs = if ih < nh
            Δs = (atm.h[ih+1] - atm.h[ih]) / cos(θ)
        else
            Δs = (atm.h[ih] - atm.h[ih-1]) / cos(θ)
        end

        # add species κ, ϵ  
        nλb = length(Iλb)
        κb  = zeros(Float64, nλb)
        ϵb  = zeros(Float64, nλb)
        for (k, val) in κbs
            @. κb += val
        end
        for (k, val) in ϵbs
            @. ϵb += val
        end

        integrate_intensity_over_Δs(Iλb, κb, ϵb, Δs, par) # 3
        push!(tt, time_ns())

        hdf5_path = if atm.h_iout[ih] == 1
            write_results_to_hdf5(par.paths, atm, ic, iθ, ih, linedata_pTNc, par.λb, Iλb, κb, ϵb, κbs, ϵbs)
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
    # loop over CO2 concentrations
    ic     = 1
    iθ, θ  = 1, 0.0
    for ic in 1:par.nbc
        # loop over angles
        for (iθ, θ) in enumerate(par.θ)

            # integrate along path
            @time cputimes = integrate_along_path(par, rdb, atm, moleculardata, linedata, ic, iθ, θ);            

            m1, m2 = size(cputimes)
            for im in 1:m2
                tim = cputimes[:,im]
                @printf("%d : %8.2e ms, %8.2e ms\n", im, sum(tim), Statistics.mean(tim))
            end 
            @printf("sum: %8.2e ms\n", sum(cputimes))
        end
    end
end

