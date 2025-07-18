using SimpleLog
using SpecialFileIO
using PhysConst
using SpecialFileIO

import PyPlot as plt
plt.pygui(true)

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
    eb = @. exp(b/λ)
    @. a / λ^5 / (exp(b/λ) - 1.0)
end

function compute_planck(T::Union{Float64,Vector{Float64}}, λ::Vector{Float64})
    IP = if isa(T, Vector{Float64})
        IPlanck = Matrix{Float64}(undef, length(λ), length(T))
        for i in eachindex(T)
            IPlanck[:,i] = planck_λ(T[i], λ)
        end
        IPlanck
    else
        planck_λ(T, λ)
    end
    IP
end

function save_planck_as_hdf5(hdf5_path::String, T::Union{Float64,Vector{Float64}}, λ::Vector{Float64}, I::Union{Vector{Float64},Matrix{Float64}})
    groups = Dict("TλI" => Dict("T" => T, "λ" => λ, "I" => I))
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end

function create_planck_spectrum(par)
    λ1 = 1.0e-6
    λ2 = 30.0*λ1
    nλ = 1000
    λP = collect(range(λ1, λ2, nλ))

    IP = planck_λ(par.surface_T, λP)
    save_planck_as_hdf5(par.paths.planck_single, par.surface_T, λP, IP)

    IPb = compute_planck(par.planck_Ts, par.λb)
    save_planck_as_hdf5(par.paths.planck_multi, par.planck_Ts, par.λb, IPb)
end

function initial_intensity(par)
    Iλb = if par.initial_intensity == :planck
        Iλb = planck_λ(par.surface_T, par.λb)
        Iλb .* (1.0 - par.albedo)
    else
        zeros(Float64, length(λλ))
    end
    save_planck_as_hdf5(joinpath(par.paths.initial_intensity), par.planck_Ts, par.λb, Iλb)
    Iλb
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