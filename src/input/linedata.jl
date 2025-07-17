using PhysConst
using SimpleLog
using Printf
using Interpolations
using DataFrames
using CSV

const TREF = 296.0 # https://hitran.org/docs/definitions-and-units/

@inline function S_T(S21r, E1, E2, β, βr, QT, QTr)
    ΔE21 = E2-E1
    S21r * QTr/QT * exp(-E1*(β+βr)) * (1.0-exp(-ΔE21*β)) / (1.0-exp(-ΔE21*βr))
end

struct LineData
    species :: Symbol
    iso   :: Vector{Int64}
    λ210  :: Vector{Float64}
    ΔE21  :: Vector{Float64}
    E1    :: Vector{Float64}
    E2    :: Vector{Float64}
    A21   :: Vector{Float64}
    B21   :: Vector{Float64}
    B12   :: Vector{Float64}
    g2    :: Vector{Float64}
    g1    :: Vector{Float64}
    S21r  :: Vector{Float64} # at TREF
    γair  :: Vector{Float64} # at TREF
    γself :: Vector{Float64} # at TREF
    nair  :: Vector{Float64}
    δair  :: Vector{Float64} # at TREF
end

function load_hitran_data(hitran_out, λmin, λmax, iso_max)
    # molec_id, local_iso_id, nu, sw, a, gamma_air, gamma_self, elower, n_air, delta_air, gp, gpp
    df0 = CSV.read(hitran_out, DataFrame)

    νmax = 1.0e-2/λmin
    νmin = 1.0e-2/λmax
    ids(x) = @. ( (x >= νmin) && (x < νmax) )
    df = df0[ids(df0.nu),:]

    iso_min = 1
    iso_max = min(maximum(df[!,:local_iso_id]), iso_max)
    isos(x) = @. ( (x >= iso_min) && (x <= iso_max) )
    df = df[isos(df.local_iso_id),:]

    sort!(df, :nu, rev=true)

    mid   = df[!,1]
    iso   = df[!,2]
    ν21   = df[!,3]                                                  # [1/cm]
    S21r  = df[!,4]                                                  # [cm^−1/(molecule ⋅ cm^−2)] = [cm / molecule], the spectral line intensity 
    A21   = df[!,5]                                                  # [[1/s]
    γair  = df[!,6]                                                  # [1 / (cm * atm)] = 1 / (1.0e-2 m * 1.01325e5 Pa)]
    γself = df[!,7]                                                  # [1/(cm * atm)]
    ν1    = df[!,8]                                                  # [1/cm]
    nair  = df[!,9]
    δair  = df[!,10]                                                 # [1/(cm * atm)]
    g2    = df[!,11]
    g1    = df[!,12]

    hc    = c_h*c_c

    ν1_m  = ν1  ./ 1.0e-2                                            # [1/m]
    ν21_m = ν21 ./ 1.0e-2                                            # [1/m]
    λ210  = 1.0 ./ ν21_m                                             # [m]
    E1    = hc  .* ν1_m                                              # [J]
    ΔE21  = hc  ./ λ210                                              # [J] 
    E2    = E1  .+ ΔE21                                              # [J]

    γair  = γair  ./ (1.0e-2 * 1.01325e5)                            # [1/(m*Pa)] , cm => m, atm => pascal
    γself = γself ./ (1.0e-2 * 1.01325e5)                            # [1/(m*Pa)] , cm => m, atm => pascal
    δair  = δair  ./ (1.0e-2 * 1.01325e5)                            # [1/(m*Pa)] , cm => m, atm => pascal

    S21r = S21r .* 1.0e-2                                            # [m]

    # Einstein coefficient of induced emission
    B21 = @. A21 * λ210^3 / (8.0*π * c_h)                            # [m^3 / s / Js] = [m^3 / J / s^2]
    # Einstein coefficient of absorption
    B12 = @. g2 / g1 * B21;

    iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair
end

function hitran_to_hdf5(species, hitran_out, hdf5, hdf5c, λmin, λmax, iso_max)
    iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair = load_hitran_data(hitran_out, λmin, λmax, iso_max)

    datasets = [
        ("iso"     , iso   ),
        ("λ210"    , λ210  ),
        ("ΔE21"    , ΔE21  ),
        ("E1"      , E1    ),
        ("E2"      , E2    ),
        ("A21"     , A21   ),
        ("B21"     , B21   ),
        ("B12"     , B12   ),
        ("g2"      , g2    ),
        ("g1"      , g1    ),
        ("S21r"    , S21r  ),
        ("γair"    , γair  ),
        ("γself"   , γself ),
        ("nair"    , nair  ),
        ("δair"    , δair  )
    ]

    save_arrays_to_hdf5(hdf5, string(species), datasets)

    data = Matrix{Float64}(undef, 14, length(λ210))
    data[ 1,:] = λ210  
    data[ 2,:] = ΔE21  
    data[ 3,:] = E1    
    data[ 4,:] = E2    
    data[ 5,:] = A21   
    data[ 6,:] = B21   
    data[ 7,:] = B12   
    data[ 8,:] = g2    
    data[ 9,:] = g1    
    data[10,:] = S21r  
    data[11,:] = γair  
    data[12,:] = γself 
    data[13,:] = nair  
    data[14,:] = δair

    save_arrays_to_hdf5(hdf5c, string(species), [("iso", iso), ("data", data)])
end

"""
    LineData(species, hitran_out, λmin, λmax, iso_max)    

    transform HITRAN data from hitran_file format (joinpath(data_dir, string(hitran_file, ".out")))
    to .npy (joinpath(input_dir, string(hitran_file, ".npy") and .hdf5 joinpath(input_dir, string(hitran_file, ".hdf5") formats, 
    using only values used in the computations

    
    species, hitran_out, λmin, λmax, iso_max = :CO2, CO2out, par.λmin, par.λmax, length(mdCO2.iso_a)
"""
function LineData(hdf5_path)
    iso_data = load_arrays_from_hdf5(hdf5_path)
    species = collect(keys(iso_data))[1]
    
    ld = if occursin("compact", hdf5_path)
        data  = iso_data[species]["data"]
        iso   = iso_data[species]["iso"]
        λ210  = data[ 1,:]
        ΔE21  = data[ 2,:]
        E1    = data[ 3,:]
        E2    = data[ 4,:]
        A21   = data[ 5,:]
        B21   = data[ 6,:]
        B12   = data[ 7,:]
        g2    = data[ 8,:]
        g1    = data[ 9,:]
        S21r  = data[10,:]
        γair  = data[11,:]
        γself = data[12,:]
        nair  = data[13,:]
        δair  = data[14,:]
        LineData(Symbol(species), iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair)
    else
        iso_data[species]["data"]
        iso   = iso_data["iso"]
        λ210  = iso_data["λ210"]                                              # m
        ΔE21  = iso_data["ΔE21"]                                              # J        
        E1    = iso_data["E1"]                                                # J
        E2    = iso_data["E2"]                                                # J
        A21   = iso_data["A21"]                                               # 1/s
        B21   = iso_data["B21"]                                               # m^3 / (J * s^2)
        B12   = iso_data["B12"]                                               # m^3 / (J * s^2)
        g2    = iso_data["g2"]                                                # 
        g1    = iso_data["g1"]                                                # 
        S21r  = iso_data["S21r"]                                              # m
        γair  = iso_data["γair"]                                              # 1 / (m * Pa)
        γself = iso_data["γself"]                                             # 1 / (m * Pa)
        nair  = iso_data["nair"]                                              # 
        δair  = iso_data["δair"]                                              # 1 / (m * Pa)
        LineData(Symbol(species), iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair)
    end
    ld
end

"""
    compute_line_emission_and_absorption_iλ(ld::LineData, Qref, Qiso, miso, c, T, N, p, iλ)
"""
function compute_line_emission_and_absorption_iλ(ld::LineData, Qref, Qiso, miso, c, T, N, p, iλ)
    dΩ = 1.0
    β  = 1.0/(c_kB * T)
    βr = 1.0/(c_kB * TREF)

    iso   = ld.iso[iλ]                                               # 
    λ210  = ld.λ210[iλ]                                              # m
    ΔE21  = ld.ΔE21[iλ]                                              # J        
    E1    = ld.E1[iλ]                                                # J
    E2    = ld.E2[iλ]                                                # J
    A21   = ld.A21[iλ]                                               # 1/s
    B21   = ld.B21[iλ]                                               # m^3 / (J * s^2)
    B12   = ld.B12[iλ]                                               # m^3 / (J * s^2)
    g2    = ld.g2[iλ]                                                # 
    g1    = ld.g1[iλ]                                                # 
    S21r  = ld.S21r[iλ]                                              # m
    γair  = ld.γair[iλ]                                              # 1 / (m * Pa)
    γself = ld.γself[iλ]                                             # 1 / (m * Pa)
    nair  = ld.nair[iλ]                                              # 
    δair  = ld.δair[iλ]                                              # 1 / (m * Pa)

    if iso > 11
        @warne mid, lid, λ210
    end

    Nspec = c*N                                                      #  [1/m^3]
    hc = c_h*c_c                                                     #  [J*m]

    # pressure shift
    λ21 = λ210 / (1.0 + δair * λ210 * p)
    Δpλ210 = λ21 - λ210

    # Lorentzian (pressure-broadened) HWHM, γ(p,T) 
    # γpT = (TREF/T)^nair * (γair * (p - pself) + γself*pself)
    γp = (TREF/T)^nair * (γair * p * (1.0 - c) + γself * p * c)      # [1/m]
    ΔλL = λ21^2 * γp                                                 # [m]

    # Doppler broadening
    ΔλG = sqrt(2.0 * c_kB * T / miso[iso]) / c_c * λ21

    # occupation numbers
    N1  = g1 * exp(- E1 * β) / Qiso[iso] * Nspec
    N2  = g2 * exp(- E2 * β) / Qiso[iso] * Nspec

    # emission [W/m^3]
    # ϵ * f(λ) * dλ * Δh                                             # [J / (m^2 * sr)]
    ϵ = hc/λ21 * A21 * dΩ/(4.0*π) *  N2                              # [J / (m^3 * sr)]

    # absorption coefficient [1]
    # κ * f(λ) * dλ * Δh                                             # [m]
    κ1 = c_h * λ21 / c_c * (N1 * B12 - N2 * B21)                     # Js * m * s/m / m^3 * m^3/(J*s^2) = 1

    S21  = S_T(S21r, E1, E2, β, βr, Qiso[iso], Qref[iso]) * Nspec    # [1/m^2]  
    κ2   = S21 * λ210^2                                              # [1]

    iso, S21, λ21, γp, ΔλL, ΔλG, N1, N2, miso[iso], ϵ, κ1, κ2
end

@doc raw"""
    compute_lines_emission_and_absorption(moleculardata, linedata, Nmolecules, T, N, p)

    Compute emission and absorption coefficients of the lines

    T - temperature
    N - atmosphere density
    p - pressure
    NCO2 - CO2 concentration
No4
"""
function compute_lines_emission_and_absorption!(ML::Matrix{Float64}, par, ld::LineData, Qref, Qiso, miso, c, T, N, p)
    Threads.@threads for iλ in eachindex(ld.λ210)
        iso, S21, λ21, γ, ΔλL, ΔλG, N1, N2, mass, ϵ, κ1, κ2 = compute_line_emission_and_absorption_iλ(ld, Qref, Qiso, miso, c, T, N, p, iλ)
        ML[:, iλ] = [iso, S21, λ21, γ, ΔλL, ΔλG, N1, N2, mass, ϵ, κ1, κ2]
    end
end
