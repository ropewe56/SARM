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
    lid   :: Vector{Int64}
    ν21   :: Vector{Float64}
    λ210  :: Vector{Float64}
    S21   :: Vector{Float64} # at TREF
    A21   :: Vector{Float64}
    γair  :: Vector{Float64} # at TREF
    γself :: Vector{Float64} # at TREF
    ΔE21  :: Vector{Float64}
    E1    :: Vector{Float64}
    E2    :: Vector{Float64}
    nair  :: Vector{Float64}
    δair  :: Vector{Float64} # at TREF
    g2    :: Vector{Float64}
    g1    :: Vector{Float64}
    B21   :: Vector{Float64}
    B12   :: Vector{Float64}
end

"""
    make_spectrum(hitran_file, data_dir, input_dir, λmin, λmax, lids, abus, masss, gis)

    transform HITRAN data from hitran_file format (joinpath(data_dir, string(hitran_file, ".out")))
    to .npy (joinpath(input_dir, string(hitran_file, ".npy") and .hdf5 joinpath(input_dir, string(hitran_file, ".hdf5") formats, 
    using only values used in the computations

"""
#outpath = paths.CO2out
#λmin, λmax = par.λmin, par.λmax
function LineData(species, outpath, λmin, λmax, iso_max)
    # molec_id, local_iso_id, nu,         sw,       a,     gamma_air, gamma_self, elower,     n_air, delta_air,  gp, gpp
    # molec_id, local_iso_id, nu        , sw       ,a    , gamma_air, gamma_self, elower    , n_air, delta_air , gp  , gpp
    # 2       , 5           , 555.003645, 7.486e-31,0.214, 0.0672   , 0.078     , 3244.5786 , 0.75 , -0.000716 , 170 , 170
    df0 = CSV.read(outpath, DataFrame)

    νmax = 1.0e-2/λmin
    νmin = 1.0e-2/λmax
    ids(x) = @.( (x >= νmin) && (x < νmax) )
    df = df0[ids(df0.nu),:]

    mid   = df[!,1]
    lid   = df[!,2]
    ν21   = df[!,3]    # 1/cm
    S21   = df[!,4]    # The spectral line intensity (cm−1/(molecule ⋅ cm−2)) cm/molecule 
    A21   = df[!,5]    # 1/s
    γair  = df[!,6]    # 1/(cm * atm)
    γself = df[!,7]    # 1/(cm * atm)
    ν1    = df[!,8]    # 1/cm
    nair  = df[!,9]
    δair  = df[!,10]   # 1/(cm * atm)
    g2    = df[!,11]
    g1    = df[!,12]

    hc    = c_h*c_c

    E1    = hc .* ν1 ./ 1.0e-2
    λ210  = 1.0e-2 ./ ν21
    ΔE21  = hc ./ λ210
    E2    = E1 .+ ΔE21
    γair  = γair  ./ (1.0e-2 * 1.0e5)  # cm => m, bar => pascal
    γself = γself ./ (1.0e-2 * 1.0e5)  # cm => m, bar => pascal
    δair  = δair  ./ (1.0e-2 * 1.0e5)  # cm => m, bar => pascal

    S21 = S21 .* 1.0e-2  # m 
    # S f(ν) = σ
    # f(ν) dν = f(λ) dλ , dν = 1/λ , f(ν) 1/λ^2 dλ = f(λ) dλ, f(ν) = λ^2 f(λ)
    # S λ^2 f(λ) = σ
    # κ = σ * N
    # σ21 = S21 * λ210^2 # [m^3]

    # Einstein coefficient of induced emission
    B21 = @. A21 * λ210^3 / (8.0*π * c_h)
    # Einstein coefficient of absorption
    B12 = @. g2 / g1 * B21;

    iso_min = 1
    iso_max = min(maximum(lid), iso_max)

    index = @.ifelse(lid >= iso_min && lid <= iso_max, true, false)
    mid   = mid[index]
    lid   = lid[index]
    ν21   = ν21[index]
    λ210  = λ210[index]
    S21   = S21[index]
    A21   = A21[index]
    γair  = γair[index]
    γself = γself[index]
    ΔE21  = ΔE21[index]
    E1    = E1[index]
    E2    = E2[index]
    nair  = nair[index]
    δair  = δair[index]
    g2    = g2[index]
    g1    = g1[index]
    B21   = B21[index]
    B12   = B12[index]


    index = sortperm(λ210)

    LineData(   species,
                lid[index],
                ν21[index],
                λ210[index],
                S21[index],
                A21[index],
                γair[index],
                γself[index],
                ΔE21[index],
                E1[index],
                E2[index],
                nair[index],
                δair[index],
                g2[index],
                g1[index],
                B21[index],
                B12[index])
end
