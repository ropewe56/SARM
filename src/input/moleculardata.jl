using PhysConst
using SimpleLog
using Printf
using Interpolations
using DataFrames
using CSV

const ch0H2O_ppm = 7966.0

struct MolecularData
    Qh       :: Vector{Vector{Float64}}
    ch       :: Vector{Float64}
    iso_id   :: Vector{Int64}
    iso_a    :: Vector{Float64}
    iso_m    :: Vector{Float64}
    gj       :: Vector{Int64}
end

function H2O_concentration(hi)
    h = reverse([84.977, 76.278, 67.577, 32.608, 41.176, 52.132, 13.792, 11.565, 8.095, 6.142, 3.77, 1.952, 0.137].*1.0e3)
    c_log10 = reverse([-5.8683, -5.5885, -5.4074, -5.3251, -5.3086, -5.2757, -5.1934, -4.6173, -3.465, -3.0864, -2.642, -2.3951, -2.0988])

    index = sortperm(h)
    h2    = h[index]
    cl = c_log10[index]
    itp = linear_interpolation(h2, cl, extrapolation_bc = Line())

    cli = itp(hi)
    ci = 10.0.^cli
    ci
end

function CO2_concentration(hi)
    h = [0.0, 10000.0, 70000.0]
    c = [1.0, 2.0/3.0, 2.0/3.0]
    itp = linear_interpolation(h, c, extrapolation_bc = Line())
    itp(hi)
end

function get_molcule_concentration(molecule, hi)
    ci = if molecule == :H2O
        H2O_concentration(hi)
    elseif molecule == :CO2
        CO2_concentration(hi)
    else
        @error molecule, "not implemented"
    end
    ci .* ppm
end


function get_concentrations(moleculardata, ih)
    ch = []    
    for im in eachindex(moleculardata)
        push!(ch, moleculardata[im].ch[ih] * ch0[im])
    end
    ch
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
function MolecularData(molecule, ch0_ppm, atm, isopath, TQmin, TQmax)
    iso_id, iso_a, iso_m, gj, qpaths = load_Isotope_file(isopath)
        
    niso = min(length(iso_id), length(qpaths))

    Qh = []
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
        push!(Qh, Q)
    end

    ch = get_molcule_concentration(molecule, hi) .* ch0_ppm
    
    MolecularData(Qh, ch, iso_id, iso_a, iso_m, gj)
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
