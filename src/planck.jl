using CommonUtils

"""
    Save wavelength and intensity to npy file

    λ - wavelengths array
    I - intensity array
    path - filepath
"""
function save_intensity_as_npy(λ::Vector{Float64}, I::Vector{Float64}, hdf5_path::String)
    n = size(I,1)
    A = zeros(Float64, 2, n)
    for i in 1:n
        A[1,i] = λ[i]
        A[2,i] = I[i]
    end
    if !(splitext(hdf5_path)[2] == ".npy")
       hdf5_path = string(hdf5_path, ".npy")
    end
    write_npy(hdf5_path, A)
end
function save_intensity_as_hdf5(hdf5_path::String, λ::Vector{Float64}, I::Vector{Float64})
    n = size(I,1)
    A = zeros(Float64, 2, n)
    for i in 1:n
        A[1,i] = λ[i]
        A[2,i] = I[i]
    end
    if !(splitext(hdf5_path)[2] == ".hdf5")
       hdf5_path = string(hdf5_path, ".hdf5")
    end
    save_array_as_hdf5(hdf5_path, A, script_dir=false)
end

function compute_and_save_planck(hdf5_path, λmin::Float64, λmax::Float64, nb_λ::Int64, T::Vector{Float64})
    λ0 = collect(LinRange(λmin, λmax, nb_λ))
    nT = length(T)
    IPlanck = Matrix{Float64}(undef, nT+1, nb_λ)
    IPlanck[1,:] = λ0
    for i in 1:nT
        IPlanck[i+1,:] = planck_λ(T[i], λ0)
    end
    if !(splitext(hdf5_path)[2] == ".hdf5")
        hdf5_path = string(hdf5_path, ".hdf5")
     end
     save_array_as_hdf5(hdf5_path, IPlanck, script_dir=false)
end

function compute_and_save_planck(hdf5_path, λmin::Float64, λmax::Float64, nb_λ::Int64, T::Float64)
    λ0 = collect(LinRange(λmin, λmax, nb_λ))
    IPlanck = planck_λ(T, λ0)
    save_intensity_as_hdf5(hdf5_path, λ0, IPlanck)
end

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

@doc raw"""
    Planck intensity at temperature temp and wavelength λ

```math
I = \dfrac{2 π  h  c^2}{λ^5} \dfrac{1}{\exp\left(\dfrac{h  c}{k_B T λ} - 1 \right)}
```
"""
@inline function planck_λ(T::Float64, λ::Vector{Float64})
    a = 2.0 * π * c_h * c_c^2
    b = c_h * c_c / (c_kB * T)
    @. a / λ^5 / (exp(b/λ) - 1.0)
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