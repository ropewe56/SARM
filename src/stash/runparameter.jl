using Parameters


const ppm = 1.0e-6
cc = [[-1.0 278.0*ppm], [-1.0 430.0*ppm], [-1.0 278.0*2.0*ppm]]
θ  = deg2rad(Vector{Float64}([0.0, 40.0, 80.0]))
planck_Ts = [288.0, 260.0, 240.0, 220.0]

@with_kw mutable struct RunParameter
    λmin                :: Float64                 = 14.0e-6
    λmax                :: Float64                 = 16.0e-6
    Δλ                  :: Float64                 = 1.0e-11
    ΔλL                 :: Float64                 = 1.0e-11
    Δλ_factor           :: Float64                 = 1.0
    surface_T           :: Float64                 = 288.0
    planck_Ts           :: Vector{Float64}         = planck_Ts
    initial_intensity   :: Symbol                  = :planck
    TQmin               :: Int64                   = 200.0
    TQmax               :: Int64                   = 300.0
    p_ref               :: Float64                 = 1.06e5
    T_ref               :: Float64                 = 288.0
    albedo              :: Float64                 = 0.3
    background          :: Float64                 = 1.0
    T_of_h              :: Bool                    = true
    N_of_h              :: Bool                    = true
    with_emission       :: Bool                    = true
    integrate           :: Bool                    = true
    θ                   :: Vector{Float64}         = θ
    c_ppm               :: Vector{Matrix{Float64}} = cc
    zmethod             :: Symbol                  = :equalnumber
    zmin                :: Float64                 = 0.0
    zmax                :: Float64                 = 70000.0
    dzmin               :: Float64                 = 10.0
    dzmax               :: Float64                 = 20000.0
    e                   :: Float64                 = 2.0
    nz                  :: Int64                   = 50
end
RunParameter()