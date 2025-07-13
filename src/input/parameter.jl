using Parameters

const ppm = 1.0e-6

θ = deg2rad.(Vector{Float64}([0.0, 40.0, 80.0]))
planck_Ts = [288.0, 260.0, 240.0, 220.0]

@with_kw mutable struct RunParameter
    λmin                :: Float64                 = 14.0e-6
    λmax                :: Float64                 = 16.0e-6
    Δλb                 :: Float64                 = 1.0e-11
    ΔλL                 :: Float64                 = 1.0e-11
    Δλ_factor           :: Float64                 = 1.0
    surface_T           :: Float64                 = 288.0
    planck_Ts           :: Vector{Float64}         = planck_Ts
    initial_intensity   :: Symbol                  = :planck
    TQmin               :: Int64                   = 200.0
    TQmax               :: Int64                   = 300.0
    κΔs_limit           :: Float64                 = 0.01
    p_ref               :: Float64                 = 1.06e5
    T_ref               :: Float64                 = 288.0
    albedo              :: Float64                 = 0.3
    background          :: Float64                 = 1.0
    T_of_h              :: Bool                    = true
    N_of_h              :: Bool                    = true
    with_emission       :: Bool                    = true
    integrate           :: Bool                    = true
    θ                   :: Vector{Float64}         = θ
    c_ppm               :: Matrix{Float64}         = zeros(Float64, 0, 0)
    hmethod             :: Symbol                  = :equalnumber
    hmin                :: Float64                 = 0.0
    hmax                :: Float64                 = 70000.0
    dhmin               :: Float64                 = 10.0
    dhmax               :: Float64                 = 20000.0
    e                   :: Float64                 = 2.0
    nh                  :: Int64                   = 50
    prealloc            :: Preallocated            = Preallocated()
end

function set_concentrations(cc)
    len = []
    for c in cc
        push!(len, length(c))
    end
    lmax = maximum(len)
    m = []
    for c in cc
        push!(m, if length(c) < lmax repeat(c, lmax) else c end)
    end
    reduce(hcat, m)' .* ppm
end

