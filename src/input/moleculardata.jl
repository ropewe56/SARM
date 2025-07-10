using PhysConst
using SimpleLog
using Printf
using Interpolations
using DataFrames
using CSV

struct MolecularData
    Qinp       
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

    CO2_Q_dir -- directory where T, Q values (HITRAN data) are
    n  -- number of T,Q pairs to make
    returns T, Q -- [description]
"""
function MolecularData(isopath, TQmin, TQmax)
    iso_id, iso_a, iso_m, gj, qpaths = load_Isotope_file(isopath)
    
    Qinp = []
    # read the partition function files
    for (i,path) in enumerate(qpaths)
        fpath = joinpath(dirname(isopath), path)        
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

        push!(Qinp, lip)
    end

    MolecularData(Qinp, iso_id, iso_a, iso_m, gj)
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
