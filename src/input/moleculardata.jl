using PhysConst
using SimpleLog
using Printf
using Interpolations
using DataFrames
using SpecialFileIO
using CSV
using JSON3

data_root() =  "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/data"

function data_files()
    DATADIR = data_root()
    d = Dict(:H2O => Dict(
                        :Q            => joinpath(DATADIR, "H2O", "H2O_Q", "H2O_Isotopes.txt"),
                        :out          => joinpath(DATADIR, "H2O", "H2O_rwfmt.csv"),
                        :hdf5         => joinpath(DATADIR, "H2O", "H2O_rwfmt.hdf5"),
                        :hdf5_compact => joinpath(DATADIR, "H2O", "H2O_rwfmt_compact.hdf5")
                ),
             :CO2 => Dict(
                        :Q            => joinpath(DATADIR, "CO2", "CO2_Q", "CO2_Isotopes.txt"),
                        :out          => joinpath(DATADIR, "CO2", "CO2_rwfmt.csv"),
                        :hdf5         => joinpath(DATADIR, "CO2", "CO2_rwfmt.hdf5"),
                        :hdf5_compact => joinpath(DATADIR, "CO2", "CO2_rwfmt_compact.hdf5")
                ),
    )
    open(joinpath(DATADIR, "data_files.json"), "w") do io
        JSON3.pretty(io, JSON3.write(d))#, ac=JSON3.AlignmentContext())
    end
    d
end

function get_data_files()
    open(joinpath(data_root(), "data_files.json"), "r") do io
        Dict(JSON3.read(io))
    end
end


struct MolecularData
    species  :: Symbol
    Qref     :: Vector{Float64} # at TREF
    Qisoh    :: Matrix{Float64}
    cnh      :: Vector{Float64}
    iso_id   :: Vector{Int64}
    iso_a    :: Vector{Float64}
    iso_m    :: Vector{Float64}
    gj       :: Vector{Int64}
end

function get_sorted_q_files(root)
    files = readdir(root)
    ii = []
    ff = []
    f = files[3]
    for f in files
        if f[1] == 'q'
            n = split(f, ".")[1][2:end]
            push!(ff, f)
            push!(ii, parse(Int64, n))
        end
    end
    index = sortperm(ii)
    paths = ff[index]
    paths
end

function load_Isotope_file(isopath)
    df = CSV.read(isopath, DataFrame, header=1)

    # mass[kg] => g/mol * mass_factor
    mass_factor = 1.0e-3/6.02214076e23

    iso_a  = df[!,:IsoAbundance]
    iso_m  = df[!,"MolarMass(g)"] .* mass_factor
    gj     = df[!,:gj]

    qpaths = get_sorted_q_files(dirname(isopath))

    iso_id = collect(range(1,length(iso_a),length(iso_a)))
    iso_id, iso_a, iso_m, gj, qpaths
end

"""
    Create lokup tables for the CO2 partition functions

    https://hitran.org/docs/iso-meta/
    global ID 	local ID 	Formula 	AFGL code 	Abundance 	        Molar Mass /g·mol-1 	Q(296 K) 	Q (full range) 	gi
    7 	        1 	        12C16O2 	626 	    0.984204 	        43.98983 	            286.09 	    q7.txt 	        1


    molec - Symbol
    atm  - Atmosphere
    isopath - isotope Q data
    TQmin, TQmax

    n  -- number of T,Q pairs to make
    returns T, Q -- [description]
"""
function MolecularData(species, atmosphere, isopath, TQmin, TQmax)
    iso_id, iso_a, iso_m, gj, qpaths = load_Isotope_file(isopath)        
    niso = min(length(iso_id), length(qpaths))
    Qhiso = []
    Qref = []
    # read the partition function files
    i = 1
    for i in 1:niso
        fpath = joinpath(dirname(isopath), qpaths[i])        
        T = []
        Q = []
        open(fpath, "r") do io
            lines = readlines(io)
            for line in lines
                spl0 = split(line, " ")
                spl = filter(x -> x != "", spl0)
                push!(T, parse(Float64, spl[1]))
                push!(Q, parse(Float64, spl[2]))
            end
        end   
        index = @. ifelse(T >= TQmin && T <= TQmax, true, false)
        T = T[index]
        Q = Q[index]
        lip = linear_interpolation(T, Q, extrapolation_bc = Line())
        Q = lip.(atmosphere.T)
        push!(Qhiso, Q)
        push!(Qref, lip(TREF))
    end
    Qisoh = reduce(hcat, Qhiso)'
    cnh = get_normalized_molecule_concentration_over_h(species, atmosphere.h)
    
    MolecularData(species, Qref, Qisoh, cnh, iso_id, iso_a, iso_m, gj)
end
mpar = par.m
spec = :CO2
function get_molecular_data(mpar, atmosphere)
    datfiles = get_data_files()
    md = Dict{Symbol,MolecularData}()
    for spec in mpar.species
        isopath, TQmin, TQmax = datfiles[spec][:Q], mpar.TQmin, mpar.TQmax
        md[spec] = MolecularData(spec, atmosphere, isopath, TQmin, TQmax)
    end
    md
end

function save_input_to_hdf5(hdf5_path, atmosphere::Atmosphere, molec_data_dict::Dict{Symbol, MolecularData})
    #h_iout:: Vector{Int64}
    #h     :: Vector{Float64}
    #p     :: Vector{Float64}
    #T     :: Vector{Float64}
    #N     :: Vector{Float64}

    #species  :: Symbol
    #Qref     :: Vector{Float64} # at TREF
    #Qisoh    :: Matrix{Float64}
    #cnh      :: Vector{Float64}
    #iso_id   :: Vector{Int64}
    #iso_a    :: Vector{Float64}
    #iso_m    :: Vector{Float64}
    #gj       :: Vector{Int64}

    datasets = Dict()
    db_path
    dbpath = joinpath(@__DIR__, "inputdata.db")
    db = SQLite.DB(dbpath)

    df_at = DataFrame(h_iout = atmosphere.h_iout, h = atmosphere.h, p = atmosphere.p, T = atmosphere.T, N = atmosphere.N)
    SQLite.load!(df2, db, "hihpTN")

    for (spec, md) in molec_data_dict
        df = 
        datasets[String(spec)] = (
            ("Qref"  , md.Qref), 
            ("Qisoh" , md.Qisoh), 
            ("cnh"   , md.cnh), 
            ("iso_id", md.iso_id), 
            ("iso_a" , md.iso_a), 
            ("iso_m" , md.iso_m), 
            ("gj"    , md.gj))
    end

    save_arrays_to_hdf5(hdf5_path, datasets; fmod="w")
    nothing
end

function load_input_from_hdf5(hdf5_path)
    groups = load_groups_as_hdf5(hdf5_path)
    atm = groups["atmosphere"]
    atmosphere = Atmosphere(atm["h_iout"], atm["h"], atm["p"], atm["T"], atm["N"])

    species = filter(x -> x != "atmosphere", collect(keys(groups)))
    molec_data = Dict()
    for spec in species
        md = groups[spec] 
        molec_data[Symbol(spec)] = MolecularData(Symbol(spec), 
            md["Qref"],
            md["Qisoh"],
            md["cnh"],
            md["iso_id"],
            md["iso_a"],
            md["iso_m"],
            md["gj"])
    end
    atmosphere, molec_data
end
