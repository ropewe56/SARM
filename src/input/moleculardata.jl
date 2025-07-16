using PhysConst
using SimpleLog
using Printf
using Interpolations
using DataFrames
using CSV


const datadir = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/data"
const H2O_Q   = joinpath(datadir, "H2O", "H2O_Q", "H2O_Isotopes.txt")
const CO2_Q   = joinpath(datadir, "CO2", "CO2_Q", "CO2_Isotopes.txt")
const H2Oout  = joinpath(datadir, "H2O", "H2O_rwfmt.out")
const CO2out  = joinpath(datadir, "CO2", "CO2_rwfmt.out")


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
    for f in files
        if occursin("q", f)
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
function MolecularData(species, atm, isopath, TQmin, TQmax)
    iso_id, iso_a, iso_m, gj, qpaths = load_Isotope_file(isopath)        
    niso = min(length(iso_id), length(qpaths))
    Qhiso = []
    Qref = []
    # read the partition function files
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
        Q = lip.(atm.T)
        push!(Qhiso, Q)
        push!(Qref, lip(TREF))
    end
    Qisoh = reduce(hcat, Qhiso)'
    cnh = get_normalized_molecule_concentration_over_h(species, atm.h)
    
    MolecularData(species, Qref, Qisoh, cnh, iso_id, iso_a, iso_m, gj)
end


function test()
    Tmin, Tmax, nT = 200.0,300.0, 100
    λmin, λmax = 14.0e-6, 16.0e-6

    datadir = "/home/wester/Projects/Julia/Climate-Energy/SARM/data"

    H2Oiso = joinpath(datadir, "H2O", "H2O_Q", "H2O_Isotopes.txt")
    CO2iso = joinpath(datadir, "CO2", "CO2_Q", "CO2_Isotopes.txt")

    H2Oout = joinpath(datadir, "H2O", "H2O_rwfmt_ISO-0-1.out")
    CO2out = joinpath(datadir, "H2O", "CO2_rwfmt_ISO-0-12_wl-12-18-mum.out")

    mdH2O = get_TQ(H2Oiso, Tmin, Tmax, nT);
    mdCO2 = get_TQ(CO2iso, Tmin, Tmax, nT);

    H2O_line_data = get_line_data(H2Oout, λmin, λmax);
    CO2_line_data = get_line_data(CO2out, λmin, λmax);
end
