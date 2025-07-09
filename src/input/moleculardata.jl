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

struct LineDataEx
    λ_ul0  
    E_l    
    E_u    
    S      
    A_ul   
    γ_a    
    γ_s    
    n_a    
    δ_a    
    g_u    
    g_l    
    iso_m  
    iso_c  
    B_ul
    B_lu
    ΔλL0
end

function get_lines(ld)
    # spectral data (HITRAN)

    @infoe @sprintf("max_isotope_id = %d", max_isotope_id)
    iso_ids = Vector{Int64}(undef,0)
    for iline in 1:nb_lines
        iso_id = spectral_data[12, iline]
        isotope_id = floor(Int64, iso_id)
        if isotope_id <= max_isotope_id
            push!(iso_ids, isotope_id)
        end
    end
    @infoe @sprintf("min_iso_id = %d, max_iso_id = %d, size(iso_ids,1) = %d", minimum(iso_ids), maximum(iso_ids), size(iso_ids,1))

    linesc = Matrix{Float64}(undef, 16, size(iso_ids,1))
    i = 1
    for iλ in 1:nb_lines
        iso_id =  spectral_data[12, iline]
        isotope_id = floor(Int64, iso_id)
        if isotope_id <= max_isotope_id

            λ_ul0  = spectral_data[ 1, iline]
            E_l    = spectral_data[ 2, iline]
            E_u    = spectral_data[ 3, iline]
            S      = spectral_data[ 4, iline]
            A_ul   = spectral_data[ 5, iline]
            γ_a    = spectral_data[ 6, iline]
            γ_s    = spectral_data[ 7, iline]
            n_a    = spectral_data[ 8, iline]
            δ_a    = spectral_data[ 9, iline]
            g_u    = spectral_data[10, iline]
            g_l    = spectral_data[11, iline]
            iso_m  = spectral_data[13, iline]
            iso_c  = spectral_data[14, iline]

            # Einstein coefficient of induced emission
            B_ul = A_ul * λ_ul0^3 / (8.0*π * c_h)
            # Einstein coefficient of absorption
            B_lu = g_u / g_l * B_ul;

            # Line pressure broadening coeffcient
            # Lorentz line width
            ΔλL0 = λ_ul0^2 * (γ_a * 1.0e5)

            linesc[i_λ_ul0, i] = λ_ul0
            linesc[i_E_l  , i] = E_l
            linesc[i_E_u  , i] = E_u
            linesc[i_S    , i] = S
            linesc[i_A_ul , i] = A_ul
            linesc[i_γ_a  , i] = γ_a
            linesc[i_γ_s  , i] = γ_s
            linesc[i_n_a  , i] = n_a
            linesc[i_δ_a  , i] = δ_a
            linesc[i_g_u  , i] = g_u
            linesc[i_g_l  , i] = g_l
            linesc[i_B_ul , i] = B_ul
            linesc[i_B_lu , i] = B_lu
            linesc[i_ΔλL0 , i] = ΔλL0
            linesc[i_iso_m, i] = iso_m
            linesc[i_iso_c, i] = iso_c

            i += 1
        end
    end

    nb_lines = size(iso_ids,1)
    linesv = Matrix{Float64}(undef, 8, nb_lines)
    LineData(iso_ids, linesc, linesv, nb_lines)
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
