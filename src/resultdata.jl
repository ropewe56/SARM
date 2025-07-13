using StaticArrays

mutable struct Results{N}
    data :: Vector{StaticVector{N, Float64}}
end
function Results(N::Int64)
    data = Vector{StaticVector{N, Float64}}(undef, 0)
    Results(data)
end

"""
    add_results(res::Results, dλ, z, NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, I_λ, vϵ, vκ, vIκ)

    add data at z

    dλ        - wavelength spacing
    z         - z
    NCO2      - CO2 concentration
    θ         - angle
    T         - temperature
    N         - gas density
    ΔλL_mean  -
    ΔλD_mean  -
    I_λ       - Vector of spectral intensity
    vϵ        - Vector of spectral emission
    vκ        - Vector of spectral absorption
    vIκ       - Vector of spectral intensity x spectral absorption
"""
function add_results!(results::Results, par::RunParameter, h, cch0, θ, T, N, ΔλL_mean, ΔλD_mean, Iλb, ϵb, κb)
    nλb = length(Iλb)

    int_ϵ1  = Vector{Float64}(undef, 2)
    int_κ1  = Vector{Float64}(undef, 2)
    int_Iκ1 = Vector{Float64}(undef, 2)
    int_ϵ2  = Vector{Float64}(undef, 2)
    int_κ2  = Vector{Float64}(undef, 2)
    int_Iκ2 = Vector{Float64}(undef, 2)
    int_ϵ3  = Vector{Float64}(undef, 2)
    int_κ3  = Vector{Float64}(undef, 2)
    int_Iκ3 = Vector{Float64}(undef, 2)

    # integrate over all wavelengths
    int_I1 = sum(Iλb) * par.Δλb

    Iκ = Vector{Vector{Float64}}(undef, 2)
    for i in eachindex(κb)
        Iκ[i] = Iλb .* κb[i]
    end

    for i in eachindex(κb)
        int_ϵ1[i]  = sum(ϵb[i]) * par.Δλb
        int_κ1[i]  = Statistics.mean(κb[i])
        int_Iκ1[i] = sum(Iκ[i]) * par.Δλb
        int_κ1[i]  = Statistics.mean(κb[i])
    end

    # integrate over nλb/6,...,nλb-nλb/6  wavelengths
    n1 = floor(Int64, nλb/6)
    n2 = nλb-n1
    for i in eachindex(κb)
        int_ϵ2[i]  = sum(ϵb[i][n1:n2]) * par.Δλb
        int_κ2[i]  = Statistics.mean(κb[i][n1:n2])
        int_Iκ2[i] = sum(Iκ[i][n1:n2]) * par.Δλb
    end
    int_I2 = sum(Iλb[n1:n2]) * par.Δλb

    # integrate over nλb/4,...,nλb-nλb/4  wavelengths
    n1 = floor(Int64, nλb/4)
    n2 = nλb-n1
    for i in eachindex(κb)
        int_ϵ2[i]  = sum(ϵb[i][n1:n2]) * par.Δλb
        int_κ2[i]  = Statistics.mean(κb[i][n1:n2])
        int_Iκ2[i] = sum(Iκ[i][n1:n2]) * par.Δλb
    end
    int_I3 = sum(Iλb[n1:n2]) * par.Δλb

    i = 1
    h = 0.0
    for i in eachindex(κb)
        push!(results.data, SA_F64[     h, θ, T, N, 
                                        cch0[i], ΔλL_mean[i], ΔλD_mean[i], 
                                        int_I1, int_ϵ1[i], int_κ1[i], int_Iκ1[i],
                                        int_I2, int_ϵ2[i], int_κ2[i], int_Iκ2[i], 
                                        int_I3, int_ϵ3[i], int_κ3[i], int_Iκ3[i]])
    end

    int_I1, int_ϵ1, int_κ1, int_Iκ1
end

"""
z, NCO2, θ, T, N, I_mean, ϵ_mean, κ_mean, ΔλL_mean, ΔλD_mean
"""
function write_result_data(res::Results, hdf5_path)
    n = size(res.data, 1)
    if n > 0
        N = length(res.data[1])
        a = zeros(Float64, N, n)
        for j in 1:n
            for i in 1:N
                a[i,j] = res.data[j][i]
            end
        end
        groups = Dict("a" => Dict("a" => a))
        save_groups_as_hdf5(hdf5_path, groups)
    else
        @warne "empty result_data!"
    end
end

#results = Vector{StaticVector{12, Float64}}(undef, 0)
#push!(results, SA_F64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
#write_result_data4("fname.npy", results)