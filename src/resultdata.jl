using Common
using StaticArrays

mutable struct Results{N}
    data :: Vector{StaticVector{N, Float64}}
end
function Results(N::Int64)
    data = Vector{StaticVector{N, Float64}}(undef, 0)
    Results(data)
end

function add_results(res::Results, dλ, z, NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, I_λ, vϵ, vκ, vIκ)
    nl = length(I_λ)
    int_I  = sum(I_λ)*dλ
    int_ϵ  = sum(vϵ)*dλ
    int_κ  = Statistics.mean(vκ)
    int_Iκ = sum(vIκ)*dλ

    n1 = floor(Int64, nl/6)
    n2 = nl-n1
    int_I1  = sum(I_λ[n1:n2])*dλ
    int_ϵ1  = sum(vϵ[n1:n2])*dλ
    int_κ1  = Statistics.mean(vκ[n1:n2])
    int_Iκ1 = sum(vIκ[n1:n2])*dλ

    n1 = floor(Int64, nl/4)
    n2 = nl-n1
    int_I2  = sum(I_λ[n1:n2])*dλ
    int_ϵ2  = sum(vϵ[n1:n2])*dλ
    int_κ2  = Statistics.mean(vκ[n1:n2])
    int_Iκ2 = sum(vIκ[n1:n2])*dλ

    push!(res.data, SA_F64[z, NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, int_I, int_ϵ, int_κ, int_Iκ, int_I1, int_ϵ1, int_κ1, int_Iκ1, int_I2, int_ϵ2, int_κ2, int_Iκ2])

    int_I, int_ϵ, int_κ, int_Iκ
end

"""
z, NCO2, θ, T, N, I_mean, ϵ_mean, κ_mean, ΔλL_mean, ΔλD_mean
"""
function write_result_data(res::Results, fname)
    n = size(res.data, 1)
    if n > 0
        N = length(res.data[1])
        a = zeros(Float64, N, n)
        for j in 1:n
            for i in 1:N
                a[i,j] = res.data[j][i]
            end
        end
        ext = splitext(fname)[2]
        if ext == ".npy"
            write_npy(fname, a)
        elseif ext == ".hdf5"
            save_array_as_hdf5(fname, a)
        else
            @warne @sprintf("File with extension [%s] cannot be not saved!", ext)
        end
    else
        @warne "empty result_data!"
    end
end

#results = Vector{StaticVector{12, Float64}}(undef, 0)
#push!(results, SA_F64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
#write_result_data4("fname.npy", results)