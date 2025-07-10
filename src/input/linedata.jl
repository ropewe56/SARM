using PhysConst
using SimpleLog
using Printf
using Interpolations
using DataFrames
using CSV

struct LineData
    mid   :: Vector{Float64}
    lid   :: Vector{Float64}
    ν21   :: Vector{Float64}
    λ210  :: Vector{Float64}
    S21   :: Vector{Float64}
    A21   :: Vector{Float64}
    γair  :: Vector{Float64}
    γself :: Vector{Float64}
    ΔE21  :: Vector{Float64}
    E1    :: Vector{Float64}
    E2    :: Vector{Float64}
    nair  :: Vector{Float64}
    δair  :: Vector{Float64}
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
function LineData(outpath, λmin, λmax)
    # molec_id, local_iso_id, nu        , sw       ,a    , gamma_air, gamma_self, elower    , n_air, delta_air , gp  , gpp
    # 2       , 5           , 555.003645, 7.486e-31,0.214, 0.0672   , 0.078     , 3244.5786 , 0.75 , -0.000716 , 170 , 170
    df0 = CSV.read(outpath, DataFrame)

    νmax = 1.0e-2/λmin
    νmin = 1.0e-2/λmax
    ids(x) = @.( (x >= νmin) && (x < νmax) )
    df = df0[ids(df0.nu),:]

    mid   = df[!,1]
    lid   = df[!,2]
    ν21   = df[!,3]
    λ210  = 1.0e-2 ./ ν21
    S21   = df[!,4]
    A21   = df[!,5]
    γair  = df[!,6] ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    γself = df[!,7] ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    ΔE21  = c_h .* c_c ./ λ210
    E1    = df[!,8] .* c_h .* c_c
    E2    = E1 .+ ΔE21
    nair  = df[!,9]
    δair  = df[!,10] ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    g2    = df[!,11]
    g1    = df[!,12]

    # Einstein coefficient of induced emission
    B21 = @. A21 * λ210^3 / (8.0*π * c_h)
    # Einstein coefficient of absorption
    B12 = @. g2 / g1 * B21;

    index = sortperm(λ210)
    LineData(   mid[index],
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
                B12)
end
