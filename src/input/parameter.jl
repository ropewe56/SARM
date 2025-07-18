using Parameters

θs = deg2rad.(Vector{Float64}([0.0, 40.0, 80.0]))
pTs = [288.0, 260.0, 240.0, 220.0]

@with_kw mutable struct RunParameter
    κΔs_limit           :: Float64                      = 0.01
    λmin                :: Float64                      = 14.0e-6
    λmax                :: Float64                      = 16.0e-6
    nλb                 :: Int64                        = 1000000
    Δλb                 :: Float64                      = 1.0e-11
    ΔλL                 :: Float64                      = 1.0e-11
    λb                  :: Vector{Float64}              = zeros(Float64, 0)
    f_Δλ_factor         :: Float64                      = 10.0
    fL_adapt            :: Symbol                       = [:none, :scale, :tail, :scaletail][4]
    fG_adapt            :: Symbol                       = [:none, :scale, :tail, :scaletail][1]
    surface_T           :: Float64                      = 288.0
    planck_Ts           :: Vector{Float64}              = pTs
    initial_intensity   :: Symbol                       = :planck
    TQmin               :: Int64                        = 200.0
    TQmax               :: Int64                        = 300.0
    albedo              :: Float64                      = 0.3
    background          :: Float64                      = 1.0
    T_of_h              :: Bool                         = true
    N_of_h              :: Bool                         = true
    with_emission       :: Bool                         = true
    integrate           :: Bool                         = true
    θ                   :: Vector{Float64}              = θs
    species             :: Vector{Symbol}               = [:H2O, :CO2]
    c_ppm               :: Dict{Symbol,Vector{Float64}} = Dict(:a => [1.0,10.0])
    nbc                 :: Int64                        = 1
    hmethod             :: Symbol                       = :equalnumber
    hmin                :: Float64                      = 0.0
    hmax                :: Float64                      = 70000.0
    dhmin               :: Float64                      = 10.0
    dhmax               :: Float64                      = 20000.0
    e                   :: Float64                      = 2.0
    nh                  :: Int64                        = 50
    prealloc            :: Preallocated                 = Preallocated()
    outdir              :: String                       = "results"
    paths               :: NamedTuple                   = NamedTuple()
end

function set_paths!(rp, outdir)
    paths = make_outpaths(outdir)
    rp.paths = paths
end

function parameter_init(par)
    par.nλb = floor(Int64, (par.λmax - par.λmin) / par.Δλb)
    par.λb  = collect(range(par.λmin, par.λmax, par.nλb))
    create_planck_spectrum(par)

    for spec in keys(par.c_ppm)
        par.c_ppm[spec][:] *= PPM 
    end
    cch0 = [par.c_ppm[k][1] for k in keys(par.c_ppm)]
    par.nbc = maximum([length(par.c_ppm[k]) for k in keys(par.c_ppm)])
end