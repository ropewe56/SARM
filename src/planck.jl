using SimpleLog
using SpecialFileIO

@doc raw"""
    Planck intensity at temperature temp and wavelength λ

```math
I = \dfrac{2 π  h  c^2}{λ^5} \dfrac{1}{\exp\left(\dfrac{h  c}{k_B T λ} - 1 \right)}
```
"""
@inline function planck(T::Float64, λ::Float64)
    a = 2.0 * π * c_h * c_c^2
    b = c_h * c_c / (c_kB * T)
    a / λ^5 / (exp(b/λ) - 1.0)
end
@inline function planck_λ(T::Float64, λ::Vector{Float64})
    a = 2.0 * π * c_h * c_c^2
    b = c_h * c_c / (c_kB * T)
    @. a / λ^5 / (exp(b/λ) - 1.0)
end

function compute_planck(T::Union{Float64,Vector{Float64}}, λ::Vector{Float64})
    IP = if isa(T, Vector{Float64})
        IPlanck = Matrix{Float64}(undef, lenght(λ), lenght(Ts))
        for i in eachinde(Ts)
            IPlanck[:,i] = planck_λ(T[i], λ)
        end
        IPlanck
    else
        planck_λ(T, λ)
    end
    IP
end

function save_planck_as_hdf5(hdf5_path::String, T::Union{Float64,Vector{Float64}}, λ::Vector{Float64}, I::Union{Vector{Float64},Matrix{Float64}})
    groups = Dict("intensity" => Dict("T" => T, "wl" => λ, "I" => I))
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end

function test_planck()
    λ = collect(LinRange(1.0e-5, 2.0e-5, 1000))
    T = collect(LinRange(200.0, 240.0, 10))
    @infoe T
    for t in T 
        f = planck_λ(t, λ)
        plot(λ*1.0e6, f)
        axis([nothing, nothing, 0.0, 1.1e7])
    end
end