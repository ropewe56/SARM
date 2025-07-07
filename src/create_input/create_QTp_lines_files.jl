using DelimitedFiles
using PhysConst
using SimpleLog
using Printf
using DataInterpolations
using DataFrames
using CSV

#using Arrow
#open(joinpath(datdir, splitext(path)[1] * ".arrow"), "w") do io
#    Arrow.write(io, df)
#end

"""
    unused
"""
function lines_to_matrix(lines)
    n = size(lines,1)
    A = []
    vmax = 0
    for i in 1:n
        sp = split(lines[i])
        for s in sp
            push!(v, parse(Float64, s))
        end
        push!(A,v)
        vmax = max(vmax, size(v, 1))
    end
    B = Matrix{Float64}(undef, vmax, n)
    for j in 1:n
        for i in 1:size(A[i],1)
            B[i,j] = A[i][j]
        end
    end
    B
end

"""
    load_out_file(datdir, path; nh=0)
    loads a HITRAN *.out file and stores the data in a *csv file
"""
function load_txt_file(datdir, path; nh=0)
    lines = open(joinpath(datdir, path), "r") do io
        split(read(io, String), "\n")[nh+1:end]
    end

    header = split(lines[1],",")
    n = length(header)
    m = length(lines[2:end])

    dat = Matrix{Float64}(undef, m, n-2)
    mid = Matrix{Int64}(undef, m, 2)
    
    for (j,line) in enumerate(lines[2:end])
        spline = split(line, ",")
        mid[j,1] = parse(Int64, spline[1])
        mid[j,2] = parse(Int64, spline[2])
        
        for (i,a) in enumerate(spline[3:end])
            if length(a) > 0
                dat[j,i] = parse(Float64, a)
            end
        end
    end

    df = DataFrame( molec_id     = mid[:,1],
                    local_iso_id = mid[:,2],                
                    nu           = dat[:,1],      
                    sw           = dat[:,2],      
                    a            = dat[:,3],     
                    gamma_air    = dat[:,4],             
                    gamma_self   = dat[:,5],              
                    elower       = dat[:,6],          
                    n_air        = dat[:,7],         
                    delta_air    = dat[:,8],             
                    gp           = dat[:,9],      
                    gpp          = dat[:,10])  

    CSV.write(joinpath(datdir, splitext(path)[1] * ".csv"), df)
end

datdir = "/home/wester/Projects/Julia/Climate-Energy/SARM/data"
path   = "H2O_rwfmt_ISO-0-1.out"
@time load_txt_file(datdir, path; nh=0)

"""
    Interpolate data onto a equidistant grid with n grid points

    x0 Vector{Float64} -- [description]
    y0 Vector{Float64} -- [description]
    n Int64 -- [description]

    (x, y) -- [description]
"""
function interpolate(x0, y0, n)
    n0 = size(x0,1)

    x = collect(LinRange(x0[1], x0[end], n))
    y = zeros(Float64, n)

    j = 1
    for i in 1:n
        xx = x[i]
        while !(x0[j] <= xx && x[i] <= x0[j+1])
            j += 1
            if (j > n0-2)
                break
            end
        end
        j = min(j, n0-2)
        v = (x[i] - x0[j]) / (x0[j+1] - x0[j])
        y[i] = (1.0 - v) * y0[j] + v * y0[j+1]
    end
    x, y
end

function load_description_file(datdir)
    lines = open(joinpath(datdir, "CO2_Q", "q7-q122-description.txt"), "r") do io
        split(read(io, String), "\n")
    end

    paths      = []
    isotope_id = []
    isotope_c  = []
    isotope_m  = []
    gis        = []

    # mass[kg] => g/mol * mass_factor
    mass_factor = 1.0e-3/6.02214076e23

    # read the description file
    for (i, line) in enumerate(lines[3:end])
        ls = split(line)
        if size(ls,1) > 5
            global_id = parse(Int64, ls[1])
            push!(isotope_id, parse(Int64, ls[2]))
            push!(isotope_c,  parse(Float64,ls[5]))
            push!(paths, ls[8])
            push!(gis, parse(Int64, ls[9]))

            mass = parse(Float64, ls[6]) * mass_factor
            push!(isotope_m, mass)
        end
    end
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

    T = []
    Q = []
    # read the partition function files
    for path in paths
        lines = readlines(joinpath(datdir, "CO2_Q", path))

        TQ = Matrix{Float64}(undef, 2, size(lines,1))
        for (i,line) in enumerate(lines)
            a = split(line)
            TQ[1, i] = parse(Float64, a[1])
            TQ[2, i] = parse(Float64, a[2])
        end

        TT = TQ[1,:]
        QQ = TQ[2,:]
        index = ifelse.(TT .< Tmax, true, false)
        push!(T, TT[index])
        push!(Q, QQ[index])
    end

    n = size(T[1],1)
    m = size(T, 1)
    TQ = zeros(Float64, m+1, n)
    TQ[1,:] = T[1]
    for i in 1:m
        TQ[i+1,:] = Q[i]
    end

    TQ, paths, isotope_id, isotope_c, isotope_m, gis
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

    lip = LinearInterpolation(p_p, h_p)
    h_p = collect(LinRange(h_p[1], h_p[end], np))
    p_p = @. lip(h_p)


    h_T = [0.0, 13.0, 17.0, 25.0, 30.0, 45.0, 50.0, 70.0] .* 1.0e3
    T_T = [288.0, 215.8, 215.7, 225.1, 233.7, 269.9, 275.7, 218.1]

    liT = LinearInterpolation(T_T, h_T)
    h_T = collect(LinRange(h_T[1], h_T[end], np))
    T_T = @. liT(h_T)

    T = zeros(Float64, 2, size(h_T,1))
    p = zeros(Float64, 2, size(h_p,1))

    T[1,:] = h_T
    T[2,:] = T_T

    p[1,:] = h_p
    p[2,:] = p_p

    h_T_path = joinpath(input_dir, "h_T.npy")
    h_p_path = joinpath(input_dir, "h_p.npy")
    write_npy(h_T_path, T)
    write_npy(h_p_path, p)

    h_T_path = joinpath(input_dir, "h_T.hdf5")
    h_p_path = joinpath(input_dir, "h_p.hdf5")
    save_array_as_hdf5(h_T_path, T, script_dir=false)
    save_array_as_hdf5(h_p_path, p, script_dir=false)
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

