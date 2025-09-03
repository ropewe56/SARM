using Parameters
using JSON3
using Dates
using Printf
using OrderedCollections

function new_subdir_path()
    d = Dates.now()
    @sprintf("%s", d)
end

function clear_subdir(result_root, subdir)
    rmpath = @sprintf("%s/*", joinpath(result_root, subdir))
    for (root, dirs, files) in walkdir(rmpath)
        for file in files
            if !occursin("", file)
                filepath = joinpath(root, file)
                rm(filepath; force=true)
            end
        end
    end
    
#    run(`bash -c "find $rmpath -type f -exec rm {} \;"`)
end

@kwdef mutable struct Paths
    outroot             :: String = ""
    intensity           :: String = ""
    spectrum            :: String = ""
    logfile             :: String = ""
    dbpath              :: String = ""
    planck_single       :: String = ""
    planck_multi        :: String = ""
    init_intensity_path :: String = ""
    input_data          :: String = ""
    molecdata_path      :: String = ""
    linedata_path       :: Dict{String,String} = Dict{String,String}()
end

function make_outpaths(result_root, subdir)
    root      = joinpath(result_root, subdir)
    intensity = joinpath(root, "intensity")
    spectrum  = joinpath(root, "spectrum")
    mkpath(root)
    mkpath(intensity)
    mkpath(spectrum)

    projroot = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl"
    molecdata_path = joinpath(projroot, "results", subdir, "input_data.hdf5")
    linedata_path  = Dict(  "CO2" => joinpath(projroot, "data/CO2/CO2_rwfmt.hdf5"), 
                            "H2O" => joinpath(projroot, "data/H2O/H2O_rwfmt.hdf5"))

    Paths(  root,
            intensity,
            spectrum,
            joinpath(root, "log.out"),
            joinpath(root, "db.sqlite3"),
            joinpath(intensity, "planck_single.hdf5"),
            joinpath(intensity, "planck_multi.hdf5"),
            joinpath(intensity, "initial_intensity.hdf5"),
            joinpath(root, "input_data.hdf5"),
            molecdata_path,
            linedata_path)
end

@kwdef mutable struct Wavelength
    λmin :: Float64 = 12.0e-6
    λmax :: Float64 = 18.0e-6
    Δλb  :: Float64 = 1.0e-11
    nλb  :: Int64   = -1
end

@kwdef mutable struct InitConditions
    initial_intensity :: Symbol = :planck  #
    planck_Ts         :: Vector{Float64} = Vector{Float64}([288.0])
    surface_T         :: Float64 = 288.0
end

@kwdef mutable struct RunParameter
    κΔs_limit        :: Float64 = 0.01
    f_Δλ_factor      :: Float64 = 10.0
    background       :: Float64 = 1.0
    θ                :: Vector{Float64} = Vector{Float64}(undef,0)
    T_of_h           :: Bool = true
    N_of_h           :: Bool = true
    f_adapt          :: Bool = true
    omit_absorb_emit :: Symbol = [:omit_none, :omit_emission, :omit_absorption][1]
end
   
@kwdef mutable struct MolecData
    TQmin    :: Float64                       = 200.0
    TQmax    :: Float64                       = 300.0
    species  :: Vector{Symbol}                = [:H2O, :CO2]
    c_ppm    :: Dict{Symbol, Vector{Float64}} = Dict(:CO2 => [1.0,10.0], :H2O => [1.0,10.0]) 
    nc       :: Int64                         = -1
end

@kwdef mutable struct Hight
    hmin    :: Float64 = 0.0
    hmax    :: Float64 = 70000.0
    dhmin   :: Float64 = 1.0
    dhmax   :: Float64 = 1000.0
    nh      :: Int64   = 500
    he      :: Float64 = 1.0
    hmethod :: Symbol  = :dh
    hout    :: Vector{Float64} = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 7000.0, 10000.0, 20000.0, 40000.0, 70000.0]
end

@kwdef mutable struct SarmParameter
    p :: Paths          = Paths()
    w :: Wavelength     = Wavelength()
    c :: InitConditions = InitConditions()
    r :: RunParameter   = RunParameter()
    m :: MolecData      = MolecData()
    h :: Hight          = Hight()
end

function rm_subdirs(result_root, subdirs)
    for sd in subdirs
        rmpath = @sprintf("%s", joinpath(result_root, sd))
        run(`rm -rf $rmpath`)
    end
end


function parameter_init_and_save(par)
    par.w.nλb = floor(Int64, (par.w.λmax - par.w.λmin) / par.w.Δλb)
    par.m.nc = maximum([length(par.m.c_ppm[k]) for k in keys(par.m.c_ppm)])

    jpath = joinpath(par.p.outroot, "parameter.json")
    to_json(jpath, par)
    JSON3.read(jpath, SarmParameter)
end

function to_json(json_path, par)
    open(json_path, "w") do io
        JSON3.pretty(io, JSON3.write(par))
    end
end

function from_json(json_path)
    jsondata = open(json_path, "r") do io
        JSON3.read(io)
    end
    dict = Dict(jsondata)
    spec = []
    for s in dict[:species]
        push!(spec, Symbol(s))
    end
    dict[:species] = spec
    cppm = Dict(dict[:c_ppm])

    cc = Dict()
    for (k,v) in cppm
        cc[Symbol(k)] = v
    end
    dict[:c_ppm] = cc

    dict[:initial_intensity] = Symbol(dict[:initial_intensity])
    dict[:hmethod] = dict[:hmethod]
    dict[:θ] = collect(dict[:θ])
    dict[:paths] = Dict(dict[:paths])

    dict[:planck_Ts] = collect(dict[:planck_Ts])

    par = RunParameter()

    for (k,v) in dict
        par.k = v
    end

    RunParameter(; dict...)
end

function make_λb(par)
    nλb = floor(Int64, (par.w.λmax - par.w.λmin) / par.w.Δλb)
    Δλb = par.w.Δλb
    λb = Vector{Float64}(undef, nλb)
    λb[1] = par.w.λmin
    for i in 2:nλb
        λb[i] = λb[i-1] + Δλb
    end
    #collect(range(par[:wavelength][:λmin], par[:wavelength][:λmax], par[:wavelength][:nλb]))
    λb
end

function show_all_par(par)
    show(IOContext(stdout, :limit=>false), MIME"text/plain"(), par)
end