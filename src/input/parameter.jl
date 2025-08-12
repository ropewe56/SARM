using Parameters
using JSON3
using Dates
using Printf
using OrderedCollections

function new_subdir_path()
    d = Dates.now()
    @sprintf("%s", d)
end

function clear_subdir(subdir)
    rmpath = @sprintf("%s/*", joinpath(OUTROOT, subdir))
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

function rm_subdirs(subdirs)
    for sd in subdirs
        rmpath = @sprintf("%s", joinpath(OUTROOT, sd))
        run(`rm -rf $rmpath`)
    end
end

function make_outpaths(subdir)
    root      = joinpath(OUTROOT, subdir)
    intensity = joinpath(root, "intensity")
    spectrum  = joinpath(root, "spectrum")
    mkpath(root)
    mkpath(intensity)
    mkpath(spectrum)

    projroot = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl"
    molecdata_path = joinpath(projroot, "results", subdir, "input_data.hdf5")
    linedata_path  = Dict(  "CO2" => joinpath(projroot, "data/CO2/CO2_rwfmt.hdf5"), 
                            "H2O" => joinpath(projroot, "data/H2O/H2O_rwfmt.hdf5"))

    OrderedDict(
        :outroot           => root,
        :intensity         => intensity,
        :spectrum          => spectrum,
        :logfile           => joinpath(root, "log.out"),
        :dbpath            => joinpath(root, "db.sqlite3"),
        :planck_single     => joinpath(intensity, "planck_single.hdf5"),
        :planck_multi      => joinpath(intensity, "planck_multi.hdf5"),
        :init_intensity_path => joinpath(intensity, "initial_intensity.hdf5"),
        :input_data        => joinpath(root, "input_data.hdf5"),
        :molecdata_path    => molecdata_path,
        :linedata_path     => linedata_path
    )
end

function get_parameter()
    OrderedDict{Symbol, Any}(
        :λmin                => 12.0e-6,
        :λmax                => 18.0e-6,
        :κΔs_limit           => 0.01,
        :Δλb                 => 1.0e-11,
        :f_Δλ_factor         => 10.0,
        :TQmin               => 200.0,
        :TQmax               => 300.0,
        :surface_T           => 288.0,
        :background          => 1.0,
        :albedo              => 0.3,
        :hmin                => 0.0,
        :hmax                => 70000.0,
        :dhmin               => 1.0,
        :dhmax               => 1000.0,
        :he                  => 1.0,
        
        :nλb                 => 1000000,
        :nh                  => 500,
        :nc                  => 2,

        :T_of_h              => true,
        :N_of_h              => true,
        :integrate           => true,

        :planck_Ts           => [288.0],
        :θ                   => [0.0],

        :species             => [:H2O, :CO2],
        :c_ppm               => Dict(:CO2 => [1.0,10.0], :H2O => [1.0,10.0]),

        :omit_absorb_emit    => [:omit_none, :omit_emission, :omit_absorption][1],
        :initial_intensity   => :planck,
        :hmethod             => :dh,
        :hout                => [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 7000.0, 10000.0, 20000.0, 40000.0, 70000.0],
        :fL_adapt            => [0, 1, 2, 3][4], # [:none, :scale, :tail, :scaletail][4],
        :fG_adapt            => [0, 1, 2, 3][1], # [:none, :scale, :tail, :scaletail][1],
    )
end

function parameter_init(par)
    par[:nλb] = floor(Int64, (par[:λmax] - par[:λmin]) / par[:Δλb])

    par[:nc] = maximum([length(par[:c_ppm][k]) for k in keys(par[:c_ppm])])

    to_json(joinpath(par[:paths][:outroot], "parameter.json"), par)
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

make_λb(par) = collect(range(par[:λmin], par[:λmax], par[:nλb]))

