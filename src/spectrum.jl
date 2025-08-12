import Statistics
using StaticArrays
using CPUTime
using Interpolations
using DataFrames

using PhysConst
using SimpleLog

function integrate_results(par, h, θ, T, N, ΔλL_mean, ΔλG_mean, Iλb, ϵb, κb)
    Δλb =  par[:Δλb]

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
    int_I = [sum(Iλb) * Δλb]

    Iκ = Dict{Symbol, Vector{Float64}}()
    for (i, spec) in enumerate(keys(κb))
        Iκ[spec] = Iλb .* κb[spec]
    end

    for (i, spec) in enumerate(keys(κb))
        int_ϵ[1][spec]  = sum(ϵb[spec]) * Δλb
        int_Iκ[1][spec] = sum(Iκ[spec]) * Δλb
        int_κ[1][spec]  = Statistics.mean(κb[spec])
    end

    # integrate over nλb/6,...,nλb-nλb/6  wavelengths
    n1 = floor(Int64, nλb/6)
    n2 = nλb-n1
    push!(int_I, sum(Iλb[n1:n2]) * Δλb)
    for (i, spec) in enumerate(keys(κb))
        int_ϵ[2][spec]  = sum(ϵb[spec][n1:n2]) * Δλb
        int_Iκ[2][spec] = sum(Iκ[spec][n1:n2]) * Δλb
        int_κ[2][spec]  = Statistics.mean(κb[spec][n1:n2])
    end

    # integrate over nλb/4,...,nλb-nλb/4  wavelengths
    n1 = floor(Int64, nλb/4)
    n2 = nλb-n1
    push!(int_I, sum(Iλb[n1:n2]) * Δλb)
    for (i, spec) in enumerate(keys(κb))
        int_ϵ[3][spec]  = sum(ϵb[spec][n1:n2]) * Δλb
        int_Iκ[3][spec] = sum(Iκ[spec][n1:n2]) * Δλb
        int_κ[3][spec]  = Statistics.mean(κb[spec][n1:n2])
    end

    int_I, int_ϵ, int_κ, int_Iκ
end

function integrated_results(λb, Iλb, ϵb, κb, ϵbs, κbs)
    Δλb   = λb[2] - λb[1]
    nλb   = length(Iλb)
    nspec = length(keys(κbs))

    n61 = floor(Int64, nλb/6)
    n62 = nλb-n61
    n41 = floor(Int64, nλb/4)
    n42 = nλb-n41

    ij = [(1,nλb), (n61, n62), (n41, n42)]
    nij = length(ij)

    int_I  = Vector{Float64}(undef, nij)
    int_ϵ  = Vector{Float64}(undef, nij)
    int_Iκ = Vector{Float64}(undef, nij)
    Iλbκb = Iλb .* κb
    for (j, (i1,i2)) in enumerate(ij)
        int_I[j]  = sum(Iλb[i1:i2]) * Δλb
        int_ϵ[j]  = sum(ϵb[i1:i2])  * Δλb
        int_Iκ[j] = sum(Iλbκb) * Δλb
    end

    Iκs = Dict{Symbol, Vector{Float64}}()
    for (i, spec) in enumerate(keys(κbs))
        Iκs[spec] = Iλb .* κbs[spec]
    end
    int_ϵs  = Dict{Symbol,Vector{Float64}}()
    mean_κs = Dict{Symbol,Vector{Float64}}()
    int_Iκs = Dict{Symbol,Vector{Float64}}()
    for spec in keys(κbs)
        int_ϵs[spec]  = Vector{Float64}(undef, nij)
        mean_κs[spec] = Vector{Float64}(undef, nij)
        int_Iκs[spec] = Vector{Float64}(undef, nij)
    end
    for spec in keys(κbs)
       for (i, (i1,i2)) in enumerate(ij)
            int_ϵs[spec][i]  = sum(ϵbs[spec][i1:i2]) * Δλb
            int_Iκs[spec][i] = sum(Iκs[spec][i1:i2]) * Δλb
            mean_κs[spec][i]  = Statistics.mean(κbs[spec][i1:i2])
        end
    end

    int_I, int_ϵ, int_Iκ, int_ϵs, mean_κs, int_Iκs
end


"""
    integrate_along_path(par, atm, moleculardata, linedata, ch0, ic, iθ, θ)

    integrate_along_path(par, rdb, atmosphere, molec_data_dict, line_data_dict, ic, iθ, θ);            

"""
function integrate_along_path(par, result_db, λb, Iλb0, atmosphere, 
                    molec_data_dict::Dict{Symbol,MolecularData}, 
                    line_data_dict::Dict{Symbol,LineData}, ic, iθ, θ)

    Iλb = copy(Iλb0)

    Δλb = par[:Δλb]
    nλb = length(λb)
    surface_T = par[:surface_T]
    T_of_h = par[:T_of_h ]
    N_of_h = par[:N_of_h ]
    nλb = length(λb)

    Tmin = par[:surface_T]
    Nmin = 1.0e30

    linedata_dict = Dict{Symbol,Matrix{Float64}}()
    for spec in par[:species]
        nλl = length(line_data_dict[spec].λ210)
        linedata_dict[spec] = Matrix{Float64}(undef, 13, nλl)
    end
    ΔλL_mean = Dict{Symbol, Float64}()
    ΔλG_mean = Dict{Symbol, Float64}()
    ϵb  = zeros(Float64, nλb)
    κb  = zeros(Float64, nλb)
    ϵbs = Dict{Symbol, Vector{Float64}}()
    κbs = Dict{Symbol, Vector{Float64}}()
    for spec in par[:species]
        ϵbs[spec] = zeros(Float64, nλb)
        κbs[spec] = zeros(Float64, nλb)
    end

    logfio = open(par[:paths][:logfile], "w")
    cputimes = []
    nh = length(atmosphere.h)

    ih = 1
    spec = :CO2

    linedata_dict = Dict{Symbol,Matrix{Float64}}()
    for spec in par[:species]
        nλl = length(line_data_dict[spec].λ210)
        linedata_dict[spec] = Matrix{Float64}(undef, 13, nλl)
    end
    ϵbs = Dict{Symbol, Vector{Float64}}()
    κbs = Dict{Symbol, Vector{Float64}}()
    ΔλL_mean = Dict{Symbol, Float64}()
    ΔλG_mean = Dict{Symbol, Float64}()
    κb  = zeros(Float64, nλb)
    ϵb  = zeros(Float64, nλb)
    for spec in par[:species]
        ϵbs[spec] = zeros(Float64, nλb)
        κbs[spec] = zeros(Float64, nλb)
    end

    for (ih,h) in enumerate(atmosphere.h)
        tt = [time_ns()]
        
        # >> 1  pressure, temperature and density at height = z
        p = atmosphere.p[ih]
        T = atmosphere.T[ih]
        N = atmosphere.N[ih]
        if T_of_h == false 
            Th = surface_T
            if N_of_h == false
                N = p / (c_kB * T)
            end
        end
        if T < Tmin Tmin = T end
        if N < Nmin Nmin = N end
        # << 1

        # >> 2
        push!(tt, time_ns())
        cihic = Dict{Symbol, Float64}()

        for spec in par[:species]
            cc = par[:c_ppm][spec]
            md   = molec_data_dict[spec]
            cihic[spec] = molec_data_dict[spec].cnh[ih] * cc[ic]

            line_data = line_data_dict[spec]
            nλl = length(line_data.λ210)
            ciso = cihic[spec]

            # 1    2          3          4          5    6    7   8    9    10  11  12 13
            # iso, miso[iso], aiso[iso], Qiso[iso], S21, λ21, γp, ΔλL, ΔλG, N1, N2, ϵ, κ 
            compute_lines_emission_and_absorption!(linedata_dict[spec], par, line_data, md.Qref, md.Qisoh[:,ih], md.iso_m, md.iso_a, ciso, T, N, p);            
        end
        # >> 2
        push!(tt, time_ns())

        # << 3
        for spec in par[:species]
            sum_over_lines!(ϵbs[spec], κbs[spec], par, λb, linedata_dict[spec])
            ΔλL_mean[spec] = Statistics.mean(linedata_dict[spec][8,:])
            ΔλG_mean[spec] = Statistics.mean(linedata_dict[spec][9,:])            
        end
        # << 3
        push!(tt, time_ns())

        # >> 4
        # step size Δs = z/cos(θ)
        Δs = if ih < nh
            Δs = (atmosphere.h[ih+1] - atmosphere.h[ih]) / cos(θ)
        else
            Δs = (atmosphere.h[ih] - atmosphere.h[ih-1]) / cos(θ)
        end

        # add species ϵ, κ
        nλb = length(Iλb)
        fill!(ϵb, 0.0)
        fill!(κb, 0.0)
        for (k, val) in ϵbs
            @. ϵb += val
        end
        for (k, val) in κbs
            @. κb += val
        end

        integrate_intensity_over_Δs(Iλb, κb, ϵb, Δs, par)
        push!(tt, time_ns())

        hdf5_path = if atmosphere.h_iout[ih] == 1
            write_results_to_hdf5(par[:paths], atmosphere, ic, iθ, ih, linedata_dict, λb, Iλb, κb, ϵb, κbs, ϵbs)
        else
            "none"
        end
        # << 4
        push!(tt, time_ns())

        # >> 5
        # add results
        int_Ij, int_ϵj, int_Iκj, int_ϵs, mean_κs, int_Iκs = integrated_results(λb, Iλb, ϵb, κb, ϵbs, κbs)        

        int_I  = sum(Iλb) * Δλb
        int_ϵ  = sum(ϵb)  * Δλb
        int_Iκ = sum(Iλb .* κb) * Δλb

        insert_into_resultdb(result_db, hdf5_path, ic, iθ, ih, atmosphere.h[ih], θ, T, N, cihic, ΔλL_mean, ΔλG_mean, 
                                        int_I, int_ϵ, int_Iκ, int_ϵs, mean_κs, int_Iκs)

        for (i, spec) in enumerate(keys(cihic))
            out = @sprintf("%s, ih = %3d, h = %12.5e, c = %12.5e, I = %12.5e, ϵ = %12.5e, Iκ = %12.5e, ΔλL = %12.5e, ΔλG = %12.5e, T = %12.5e, N = %12.5e",
                                spec, ih, atmosphere.h[ih], cihic[spec], int_I, int_ϵ, int_Iκ, ΔλL_mean[spec], ΔλG_mean[spec], T, N)
            @infoe out
        end        
        # << 5
        push!(tt, time_ns())

        dt = tt[2:end] - tt[1:end-1]
        push!(cputimes, [Float64(x).*1.0e-6 for x in dt])
    end  # lop over z ih
    
    CPUt = reduce(hcat, cputimes)' .* 1.0e-3
    m1, m2 = size(CPUt)
    im = 1
    for im in 1:m2
        tim = CPUt[:,im]
        @printf("%d : sum = %8.2e s\n", im, sum(tim))
        #tmin, tmax = extrema(tim)
        #@printf("%d : sum = %8.2e s, mean = %8.2e s, min = %8.2e s, max = %8.2e s\n", im, sum(tim), Statistics.mean(tim), tmin, tmax)
    end 
    @printf("totalsum= %8.2e s\n", sum(CPUt))

end

"""
    function integrate(par, atm::Atmosphere, moleculardata::Vector{MolecularData}, linedata::Vector{LineData})
No1
"""
function integrate(par, result_db, λb, Iλb, atmosphere::Atmosphere, molec_data_dict::Dict{Symbol,MolecularData},  
                        line_data_dict::Dict{Symbol,LineData})
    ic     = 1
    iθ, θ  = 1, 0.0
    for ic in 1:par[:nc], (iθ, θ) in enumerate(par[:θ])
        @time integrate_along_path(par, result_db, λb, Iλb, atmosphere, molec_data_dict, line_data_dict, ic, iθ, θ);            
    end
end

