#using DelimitedFiles
using PhysConst
using SimpleLog
using Printf
using DataInterpolations
using DataFrames
using CSV


datadir = "/home/wester/Projects/Julia/Climate-Energy/SARM/data"
opath = "H2O_rwfmt_ISO-0-1.out"
qpath = joinpath("CO2_Q", "CO2_q7-q122.txt")

CSV.read(joinpath(datadir,opath), DataFrame)
CSV.read(joinpath(datadir,qpath), DataFrame, header=1)

lines = open(joinpath(datdir, "CO2_Q", "q7-q122-description.txt"), "r") do io

"""
    Interpolate data onto a equidistant grid with n grid points

    x0 Vector{Float64} -- [description]
    y0 Vector{Float64} -- [description]
    n Int64 -- [description]

    (x, y) -- [description]
"""
function interpolate(x0, y0, n)
    lp = LinearInterpolation(x0, y0)
    x = collect(LinRange(x0[1], x0[end], n))
    x, lp.(x)
end

function load_description_file(datadir)
    qpath = joinpath("CO2_Q", "CO2_q7-q122.txt")
    df = CSV.read(joinpath(datadir,qpath), DataFrame, header=1)

    # mass[kg] => g/mol * mass_factor
    mass_factor = 1.0e-3/6.02214076e23

    isotope_id = df[!,:localID]
    isotope_c  = df[!,:Abundance]
    isotope_m  = df[!,"MolarMass/g·mol-1"] .* mass_factor
    paths      = df[!,"Q(fullrange)"]
    gis        = df[!,:gi]

    paths, isotope_id, isotope_c, isotope_m, gis       
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
function make_lookup_for_Q(datdir; Tmax=300.0)
    paths, isotope_id, isotope_c, isotope_m, gis = load_description_file(datdir)

    df = CSV.read(joinpath(datadir, "CO2_Q", paths[1]), DataFrame, headers=0)

    # read the partition function files
    for path in paths[2:end]
        df2 = CSV.read(joinpath(datadir, "CO2_Q", path), DataFrame, headers=0)
        append!(df, df2)
    end
    T = df[!,1]
    Q = df[!,2]
    index = @. ifelse(T >= Tmin && T <= Tmax, true, false)
    T = T[index]
    Q = Q[index]

    n = size(T[1],1)
    m = size(T, 1)
    TQ = zeros(Float64, m+1, n)
    TQ[1,:] = T[1]
    for i in 1:m
        TQ[i+1,:] = Q[i]
    end

    TQ, paths, isotope_id, isotope_c, isotope_m, gis
end

struct HPT
    hT :: Vector{Float64} 
    T  :: Vector{Float64}
    hp :: Vector{Float64} 
    p  :: Vector{Float64}
end
"""
    make_T_p_over_height(input_dir; np = 100)

    height dependent digitized values of T and p http://climatemodels.uchicago.edu/modtran/
    create linear interpolation objects and interpolate values onto an equidistamt grid of n points

    save interpolate values:
        joinpath(input_dir, "h_T.npy")
        joinpath(input_dir, "h_p.npy")
        joinpath(input_dir, "h_T.hdf5")
        joinpath(input_dir, "h_p.hdf5")

"""
function make_T_p_over_height(input_dir; np = 100)
    h1 = [   0.0,   1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0,  11.0,  12.0,
            13.0,  14.0,  15.0,  16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 70.0]
    p1 = [1013.0, 902.0, 802.0, 710.0, 628.0, 554.0, 487.0, 426.0, 372.0, 324.0, 281.0, 243.0, 209.0,
           179.0, 153.0, 130.0, 111.0, 95.0, 81.2, 69.5, 59.5, 51.0, 43.4, 27.7, 13.2, 6.52, 3.33, 1.76, 0.951, 0.067]

    h_p  = h1 .* 1.0e3
    p_p  = p1 .* 1.0e2
    h_p, p_p = interpolate(p_p, h_p, np)


    h_T = [0.0, 13.0, 17.0, 25.0, 30.0, 45.0, 50.0, 70.0] .* 1.0e3
    T_T = [288.0, 215.8, 215.7, 225.1, 233.7, 269.9, 275.7, 218.1]
    h_T, T_T = interpolate(h_T, T_T, np)

    hpt = HPT(h_T, T_T, h_p, p_p)

    hpt
end

"""
    make_spectrum(hitran_file, data_dir, input_dir, λmin, λmax, lids, abus, masss, gis)

    transform HITRAN data from hitran_file format (joinpath(data_dir, string(hitran_file, ".out")))
    to .npy (joinpath(input_dir, string(hitran_file, ".npy") and .hdf5 joinpath(input_dir, string(hitran_file, ".hdf5") formats, 
    using only values used in the computations

"""
function make_spectrum(hitran_file, data_dir, input_dir, λmin, λmax, lids, abus, masss, gis)
    hpath = joinpath(data_dir, string(hitran_file, ".out"))
    (m, ms) = readdlm(hpath, ',', '\n'; header=true)

    mid = m[:,1]
    iid = m[:,2]
    ν   = m[:,3]
    ν   = ν ./ 1.0e-2
    λ   = vec(1.0 ./ ν)

    index = vec((λ .>= λmin) .& (λ .< λmax))

    λ   = λ[index]
    iid = iid[index]

    S   = m[index,  4]
    A   = m[index,  5]
    γ_a = m[index,  6]
    γ_s = m[index,  7]
    ν_l = m[index,  8]
    n_a = m[index,  9]
    δ_a = m[index, 10]
    g_u = m[index, 11]
    g_l = m[index, 12]

    ν_l    = ν_l ./ 1.0e-2           # cm => m, bar => pascal
    γ_a    = γ_a ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    γ_s    = γ_s ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    δ_a    = δ_a ./ 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    ΔE_ul  = c_h .* c_c ./ λ
    E_l    = c_h .* c_c .* ν_l
    E_u    = E_l .+ ΔE_ul

    n = size(A,1)
    c = zeros(Float64, 14, n)

    c[  1, :] = λ
    c[  2, :] = E_l
    c[  3, :] = E_u
    c[  4, :] = S
    c[  5, :] = A
    c[  6, :] = γ_a
    c[  7, :] = γ_s
    c[  8, :] = n_a
    c[  9, :] = δ_a
    c[ 10, :] = g_u
    c[ 11, :] = g_l

    i = argmax(c[4,:])
    # 0 1 2 3 4 5 6 7 8 9 10 11
    # 1 2 3 4 5 6 7 8 9 0 11 12

    itoj = [9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 10, 11] .+ 1

    for i in 1:n
        ii = Int64(iid[i])+1
        j = itoj[ii]
        c[12,i] = j
        c[13,i] = masss[j]
        c[14,i] = abus[j]
    end

    write_npy(joinpath(input_dir, string(hitran_file, ".npy")), c)
    save_array_as_hdf5(joinpath(input_dir, string(hitran_file, ".hdf5")), c,script_dir=false)

    for j in 30:50
        s = []
        for i in 1:14
            push!(s, @sprintf("%12.6e", c[i, j]))
        end
    end
end

"""
    create_data_files(data_dir, input_dir)

    create files containing:
         CO2 partition functions
         spectral data
"""
function create_data_files(data_dir, input_dir)
    # h_T_path, npy, hdf5
    # h_p_path, npy, hdf5
    make_T_p_over_height(input_dir)

    CO2_Q_dir = joinpath(data_dir, "CO2_Q")
    TQ, paths, isotope_id, isotope_c, isotope_m, gis = make_lookup_for_Q(CO2_Q_dir, Tmax = 300.0)

    TQ_path = joinpath(input_dir, "T_Q.npy")
    write_npy(TQ_path, TQ)
    TQ_path = joinpath(input_dir, "T_Q.hdf5")
    save_array_as_hdf5(TQ_path, TQ, script_dir=false)

    λmin, λmax = 1.19e-5, 1.81e-5
    hitran_file = "CO2_rwfmt_ISO-0-12_wl-12-18-mum"
    # CO2_rwfmt_ISO-0-12_wl-12-18-mum, npy, hdf5
    make_spectrum(hitran_file, data_dir, input_dir, λmin, λmax, isotope_id, isotope_c, isotope_m, gis)
end

