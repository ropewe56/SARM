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
    mid   :: Vector{Float64}
    lid   :: Vector{Float64}
    ν21   :: Vector{Float64}
    λ210  :: Vector{Float64}
    S21   :: Vector{Float64}
    A21   :: Vector{Float64}
    γair  :: Vector{Float64}
    γself :: Vector{Float64}
    ΔE21  :: Vector{Float64}
    E1    :: Vector{Float64}
    E2    :: Vector{Float64}
    nair  :: Vector{Float64}
    δair  :: Vector{Float64}
    g2    :: Vector{Float64}
    g1    :: Vector{Float64}
    B21   :: Vector{Float64}
    B12   :: Vector{Float64}
end

function get_melcule_and_isotope_ids()
    @infoe @sprintf("max_isotope_id = %d", max_isotope_id)
    iso_ids = Vector{Int64}(undef,0)
    for iline in 1:nb_lines
        iso_id = spectral_data[12, iline]
        isotope_id = floor(Int64, iso_id)
        if isotope_id <= max_isotope_id
            push!(iso_ids, isotope_id)
        end
    end
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
function get_line_data(outpath, λmin, λmax)
    # molec_id, local_iso_id, nu        , sw       ,a    , gamma_air, gamma_self, elower    , n_air, delta_air , gp  , gpp
    # 2       , 5           , 555.003645, 7.486e-31,0.214, 0.0672   , 0.078     , 3244.5786 , 0.75 , -0.000716 , 170 , 170
    df0 = CSV.read(outpath, DataFrame)

    νmax = 1.0e-2/λmin
    νmin = 1.0e-2/λmax
    ids(x) = @.( (x >= νmin) && (x < νmax) )
    df = df0[ids(df0.nu),:]

    mid   = df[!,1]
    lid   = df[!,2]
    ν21   = df[!,3]
    λ210  = 1.0e-2 ./ ν
    S21   = df[!,4]
    A21   = df[!,5]
    γair  = df[!,6] ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    γself = df[!,7] ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    ΔE21  = c_h .* c_c ./ λ210
    E1    = df[!,8] .* c_h .* c_c
    E2    = E1 .+ ΔE21
    nair  = df[!,9]
    δair  = df[!,10] ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    g2    = df[!,11]
    g1    = df[!,12]

    # Einstein coefficient of induced emission
    B21 = A21 * λ210^3 / (8.0*π * c_h)
    # Einstein coefficient of absorption
    B12 = g2 / g1 * B21;

    LineData(   mid  ,
                lid  ,
                ν21  ,
                λ210 ,
                S21  ,
                A21  ,
                γair ,
                γself,
                ΔE21 ,
                E1   ,
                E2   ,
                nair ,
                δair ,
                g2   ,
                g1   ,
                B21  ,
                B12)
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
