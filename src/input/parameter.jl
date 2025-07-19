using Parameters
using JSON3
using Dates
using Printf

const OUTROOT = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results"

function make_outpaths()
    d = Dates.now()
    subdir = @sprintf("%s", d)[1:16]

    root      = joinpath(OUTROOT, subdir)
    intensity = joinpath(root, "intensity")
    spectrum  = joinpath(root, "spectrum")
    atm       = joinpath(root, "atm")
    mkpath(root)
    mkpath(intensity)
    mkpath(spectrum)
    mkpath(atm)

    Dict(
        :outroot           => root,
        :intensity         => intensity,
        :spectrum          => spectrum,
        :atm               => atm,
        :logfile           => joinpath(root, "log.out"),
        :dbpath            => joinpath(root, "db.sqlite3"),
        :planck_single     => joinpath(intensity, "planck_single.hdf5"),
        :planck_multi      => joinpath(intensity, "planck_multi.hdf5"),
        :initial_intensity => joinpath(intensity, "initial_intensity.hdf5")
    )
end


function get_parameter()
    Dict{Symbol, Any}(
        :κΔs_limit           => 0.01,
        :λmin                => 14.0e-6,
        :λmax                => 16.0e-6,
        :nλb                 => 1000000,
        :Δλb                 => 1.0e-11,
        :ΔλL                 => 1.0e-11,
        :λb                  => zeros(Float64, 0),
        :f_Δλ_factor         => 10.0,
        :fL_adapt            => [:none, :scale, :tail, :scaletail][4],
        :fG_adapt            => [:none, :scale, :tail, :scaletail][1],
        :surface_T           => 288.0,
        :planck_Ts           => [288.0],
        :initial_intensity   => :planck,
        :TQmin               => 200.0,
        :TQmax               => 300.0,
        :albedo              => 0.3,
        :background          => 1.0,
        :T_of_h              => true,
        :N_of_h              => true,
        :with_emission       => true,
        :integrate           => true,
        :θ                   => [0.0],
        :species             => [:H2O, :CO2],
        :c_ppm               => Dict(:a => [1.0,10.0]),
        :nbc                 => 1,
        :hmethod             => :equalnumber,
        :hmin                => 0.0,
        :hmax                => 70000.0,
        :dhmin               => 10.0,
        :dhmax               => 20000.0,
        :e                   => 2.0,
        :nh                  => 50,
        :outdir              => "results",
        :paths               => Dict()
    )
end

function parameter_init(par)
    par[:nλb] = floor(Int64, (par[:λmax] - par[:λmin]) / par[:Δλb])
    λb  = collect(range(par[:λmin], par[:λmax], par[:nλb]))
    create_planck_spectrum(par, λb)

    for spec in keys(par[:c_ppm])
        par[:c_ppm][spec][:] *= PPM 
    end
    cch0 = [par[:c_ppm][k][1] for k in keys(par[:c_ppm])]
    par[:nbc] = maximum([length(par[:c_ppm][k]) for k in keys(par[:c_ppm])])

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

    dict[:fL_adapt] = Symbol(dict[:fL_adapt])
    dict[:fG_adapt] = Symbol(dict[:fG_adapt])
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

#json_path = "jsonpath.json"
#to_json("jsonpath.json", par)
#par = from_json("jsonpath.json")

using Parameters
@with_kw mutable struct Dog
    name::String
    breed::Symbol = :husky
end

d = Dict([:name => "x", :breed => :wolf])
Dog(; d...)