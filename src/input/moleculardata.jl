using PhysConst
using SimpleLog
using Printf
using Interpolations
using DataFrames
using CSV

struct MolecularData
    Tq       :: Vector{Float64}
    Qs       :: Vector{Vector{Float64}}
    iso_id   :: Vector{Int64}
    iso_a    :: Vector{Float64}
    iso_m    :: Vector{Float64}
    gj       :: Vector{Int64}
end

struct LineData
    mid:: Vector{Int64}
    lid:: Vector{Int64}
    λ  :: Vector{Float64}
    El :: Vector{Float64}
    Eu :: Vector{Float64}
    S  :: Vector{Float64}
    A  :: Vector{Float64}
    γa :: Vector{Float64}
    γs :: Vector{Float64}
    na :: Vector{Float64}
    δa :: Vector{Float64}
    gu :: Vector{Float64}
    gl :: Vector{Float64}
end

function load_description_file(datadir, subdir, qname)
    qpath = joinpath(datadir, subdir, qname)
    df = CSV.read(qpath, DataFrame, header=1)

    # mass[kg] => g/mol * mass_factor
    mass_factor = 1.0e-3/6.02214076e23

    isotope_id = df[!,:localID]
    isotope_a  = df[!,:Abundance]
    isotope_m  = df[!,"MolarMass/g·mol-1"] .* mass_factor
    qpaths     = df[!,"Q(fullrange)"]
    gis        = df[!,:gi]

    isotope_id, isotope_a, isotope_m, gis, qpaths
end

function load_Isotope_file(iso)
    df = CSV.read(iso, DataFrame, header=1)

    # mass[kg] => g/mol * mass_factor
    mass_factor = 1.0e-3/6.02214076e23

    iso_id = df[!,:Molecule]
    iso_a  = df[!,:IsoAbundance]
    iso_m  = df[!,"MolarMass(g)"] .* mass_factor
    gj     = df[!,:gj]

    files = readdir(joinpath(datadir, subdir))
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

    iso_id, iso_a, iso_m, gj, paths
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
function get_TQ(iso, Tmin, Tmax, nT)
    #isotope_id, isotope_c, isotope_m, gis, qpaths = load_description_file(datadir, subdir, qname)

    iso_id, iso_a, iso_m, gj, qpaths = load_Isotope_file(iso)
    
    Qs = []
    Tq = collect(range(Tmin, Tmax, nT))
    # read the partition function files
    for path in qpaths
        fpath = joinpath(datadir, subdir, path)
        
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
        index = @. ifelse(T >= Tmin && T <= Tmax, true, false)
        T = T[index]
        Q = Q[index]

        lip = linear_interpolation(T, Q, extrapolation_bc = Line())
        Q = lip.(Tq)
        push!(Qs, Q)
    end

    MolecularData(Tq, Qs, iso_id, iso_a, iso_m, gj)
end

"""
    make_spectrum(hitran_file, data_dir, input_dir, λmin, λmax, lids, abus, masss, gis)

    transform HITRAN data from hitran_file format (joinpath(data_dir, string(hitran_file, ".out")))
    to .npy (joinpath(input_dir, string(hitran_file, ".npy") and .hdf5 joinpath(input_dir, string(hitran_file, ".hdf5") formats, 
    using only values used in the computations

"""
function get_line_data(Mout, λmin, λmax)
    df0 = CSV.read(Mout, DataFrame)

    νmax = 1.0e-2/λmin
    νmin = 1.0e-2/λmax
    ids(x) = @.( (x >= νmin) && (x < νmax) )
    df = df0[ids(df0.nu),:]

    mid = df[!,1]
    lid = df[!,2]
    ν   = df[!,3]
    λ   = 1.0e-2 ./ ν

    S   = df[!,  4]
    A   = df[!,  5]
    γ_a = df[!,  6]
    γ_s = df[!,  7]
    ν_l = df[!,  8]
    n_a = df[!,  9]
    δ_a = df[!, 10]
    g_u = df[!, 11]
    g_l = df[!, 12]

    ν_l   = ν_l ./ 1.0e-2           # cm => m, bar => pascal
    γ_a   = γ_a ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    γ_s   = γ_s ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    δ_a   = δ_a ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    ΔE_ul = c_h .* c_c ./ λ
    E_l   = c_h .* c_c .* ν_l
    E_u   = E_l .+ ΔE_ul

    LineData(mid, lid, λ, E_l, E_u, S, A, γ_a, γ_s, n_a, δ_a, g_u, g_l)
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
