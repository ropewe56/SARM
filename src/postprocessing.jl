using QuadGK
using JSON
using PyPlot
pygui(false)
pygui(:qt5)

using Common
using PhysConst

"""

    compute_intensity(z, I, κ, ϵ)

    z -
    I -
    κ -
    ϵ -

    returns: intensity
"""
function compute_intensity(z, I, κ, ϵ)
    dz = 100.0
    exp_κ = exp(-κ * dz)
    a     = ϵ / κ * (1.0 - exp_κ)
    b     = ϵ * dz
    I_    = zeros(Float64, size(z,1))
    I_[1] = I[1]
    limit = 0.01
    ee   = 0.0
    for i in 1:size(z,1)-1
        ee = ee + ϵ[i] * dz
        if κ[i] * dz < limit
            I_[i+1] = I_[i] * exp_κ[i] + ϵ[i] * dz
        else
            I_[i+1] = I_[i] * exp_κ[i] + a[i]
        end
    end
    return I_
end


"""
    Read all result data from the log.log files
"""
mutable struct DataSets
    dirs      :: Vector{String}
    datasets  :: Dict
    parameter :: Dict
end

function makeDataSets(out_root::String, subdir::String, subsubdirs::Vector{String})
    dirs = []
    for ssd in subsubdirs
        push!(dirs, joinpath(out_root, subdir, ssd))
    end
    sort!(dirs)

    datasets = Dict()
    parameter = Dict()
    for (id, dd) in enumerate(dirs)
        datasets[id], parameter[id] = load_results(dd)
    end
    DataSets(dirs, datasets, parameter)
end

get_size(d::DataSets) = size(d,1)

"""
# iN = 1, iθ = 1, NCO2 =  4.00000e+02, θ =  0.00000e+00

iz =   0, z =  0.00000e+00,  NCO2 =  4.00000e+02, θ =  0.00000e+00, I =  1.10178e+02, T =  2.88000e+02, N =  2.54761e+25, ϵ =  0.00000e+00, κ =  0.00000e+00, ΔλL =  0.00000e+00, ΔλD =  0.00000e+00
"""
function decode_line(line)
    l1 = split(line, ",")
    keys   = Vector{String}(undef,0)
    values = Vector{Float64}(undef,0)
    for (i, a) in enumerate(l1)
        b = split(a, "=")
        push!(keys, b[1])
        push!(values, parse(Float64, b[2]))
    end
    keys, values
end

function load_results(out_root)
    files_and_dirs = readdir(out_root)
    data = []
    for f in files_and_dirs
        if occursin("result_", f)
            fs = split(splitext(f)[1], "_")
            iN = parse(Int64, fs[2])
            iθ = parse(Int64, fs[3])
            dat = load_array_as_hdf5(joinpath(out_root, f), group="A", dataset="A")
            @infoe @sprintf("%d  %d  %s", iN, iθ, f)
            push!(data, (iN, iθ, dat))
        end
    end

    parameter = JSON.parsefile(joinpath(out_root, "input.json"))#; dicttype=Dict, inttype=Int64, use_mmap=true)
    data, parameter
end

"""[summary]

Arguments:
    id {[type]} -- [description]
    out_root {[type]} -- [description]
"""
function read_data(out_root)
    parameter = JSON.parsefile(joinpath(out_root, "input.json"))#; dicttype=Dict, inttype=Int64, use_mmap=true)

    lines = readlines(joinpath(out_root, "log.log"), keep=false) # dont keeep \n

    ls = split(lines[1])

    # 2  3  700  10
    # # iN = 1, iθ = 1, NCO2 =  4.00000e+02, θ =  0.00000e+00

    nN   = parse(Int64, ls[1])     # nb NCO2
    nθ   = parse(Int64, ls[2])     # nb theta
    nz   = parse(Int64, ls[3])+1   # nb z, there are two iz = 0 rows in log.log file
    ncol = parse(Int64, ls[4])     # nb columns
    #        0     1  2  3  4  5  6  7  8    9
    #    iz, z, NCO2, θ, I, T, N, ϵ, κ, ΔλL, ΔλD
    #iz =   0, z =  0.00000e+00,  NCO2 =  4.00000e+02, θ =  0.00000e+00, I =  1.10178e+02, T =  2.88000e+02, N =  2.54761e+25, ϵ =  0.00000e+00, κ =  0.00000e+00, ΔλL =  0.00000e+00, ΔλD =  0.00000e+00

    # 3  2  701 10
    M = zeros(Float64, nN, nθ, nz, ncol+1)
    iN = 1
    iθ = 1
    inext = 0
    for (i,line) in enumerate(lines[2:end])
        #println(line, "  ", line[1] == "#" )
        if line[1] == '#' # indicates next run
            k,v = decode_line(line[2:end])
            iN = floor(Int64, v[1])
            iθ = floor(Int64, v[2])
            inext = i+1
        else
            k,v = decode_line(line)
            iz = floor(Int64, v[1])
            M[iN, iθ, iz, :] = v[2:end]
            if i == inext
                @infoe @sprintf("%d  %d  %d  %s", iN, iθ, iz, v[2:end])
            end
        end
    end
    M, parameter
end

"""
    Get data

    id {int} -- subsubdir
    iN {int} -- NCO2
    iθ {int} -- theta
    ie {int} -- column
        0     1  2  3  4  5  6  7  8    9
    iz, z, NCO2, θ, I, T, N, ϵ, κ, ΔλL, ΔλD

    array -- a column
"""
function get(d::DataSets, id, iN, iθ, ie)
    ds = d.datasets[id]
    @infoe @sprintf("%d  %d  %d  %d  %d  %d  %d  %d", id, iN, iθ, ie, ds[1][1], ds[1][2], ds[2][1], ds[2][2])
    println(ds[1][1], "  ", ds[1][2], "  ", ds[2][1], "  ", ds[2][2])
    #ds = Vector{Tuple{Int64, Int64, Matrix}}
    #       NCO2, theta, z, icol
    i = -1
    for iv in eachindex(ds)
        if ds[iv][1] == iN && ds[iv][2] == iθ
            i = iv
        end
    end
    if i > -1
        #@infoe @sprintf("%d  %d  %d  %d  %s", id, iN, iθ, ie, ds[i][3][ie,1])
        return ds[i][3][ie,:]
    end
    @warne ds[i]
end

# indices
const ciz  = 1
const ciC  = 2
const ciθ  = 3
const ciT  = 4
const ciN  = 5
const ciI  = 8
const ciϵ  = 9
const ciκ  = 10
const ciIκ = 11
const ci1I  = 12
const ci1ϵ  = 13
const ci1κ  = 14
const ci1Iκ = 15
const ci2I  = 16
const ci2ϵ  = 17
const ci2κ  = 18
const ci2Iκ = 19
# 1  2     3  4  5  6         7         8      9      10     11      12      13      14      15       16      17      18      19
# z, NCO2, θ, T, N, ΔλL_mean, ΔλD_mean, int_I, int_ϵ, int_κ, int_Iκ, int_I1, int_ϵ1, int_κ1, int_Iκ1, int_I2, int_ϵ2, int_κ2, int_Iκ2]
"""
    Integrate intensity ove angle theta

    data_sets {[type]} -- [description]
    id {[type]} -- [description]
    nθ {[type]} -- [description]
    iN {[type]} -- [description]
    NCO2 {[type]} -- [description]
    color1 {[type]} -- [description]
    color2 {[type]} -- [description]

    [type] -- [description]
"""
function integrate_intensity(data_sets, id, nθ, iN, iθ, iI, NCO2, color1, color2)
    θ_0 = get(data_sets, id, iN, 1, iθ)[1]   # theta
    θ_1 = get(data_sets, id, iN, 2, iθ)[1]   # theta
    θ_2 = get(data_sets, id, iN, 3, iθ)[1]   # theta
    I_0 = get(data_sets, id, iN, 1, iI)[end] # intensity at TOA
    I_1 = get(data_sets, id, iN, 2, iI)[end] # intensity at TOA
    I_2 = get(data_sets, id, iN, 3, iI)[end] # intensity at TOA

    #@infoe @sprintf("%s %s  %s", I_0, I_1, I_2)
    # qubic approximation of I(θ)
    R1 = I_1 - I_0
    R2 = I_2 - I_0

    a0 = I_0
    det = θ_1^2 * θ_2^3 - θ_1^3 * θ_2^2
    a2 = (R1 * θ_2^3 - R2 * θ_1^3) / det
    a3 = (R2 * θ_1^2 - R1 * θ_2^2) / det

    #x, w = gausslegendre(12)
    c1, err = quadgk(x -> cos(x)*sin(x)    , 0, π*0.5, rtol=1e-8)
    c2, err = quadgk(x -> cos(x)*sin(x)*x^2, 0, π*0.5, rtol=1e-8)
    c3, err = quadgk(x -> cos(x)*sin(x)*x^3, 0, π*0.5, rtol=1e-8)

    θ = collect(LinRange(0.0, 0.5*π, 100))
    I = @. a0 + a2 * θ^2 + a3 * θ^3

    # plot
    #plot(θ, I, color1, label='%d ppm, cubic approximation' % NCO2)
    #plot([θ_0, θ_1, θ_2], [I_0, I_1, I_2], color2+'o', label="%d ppm, computed" % NCO2)
    #xlabel("angle θ [rad]")
    #ylabel("TOA flux I(θ) [W/m²]")
    #legend(loc='best')

    # integrated intensity
    Iint = 2.0*π * (a0*c1[1] + a2*c2[1] + a3*c3[1])
    @infoe Iint
    return Iint
end

function plot_results(data_sets, id1, id2, idl)
    with_emission = data_sets.parameter[id]["run_parameter"]["with_emission"]

    plotfunc = plot
    tickfunc = xticks

    root = joinpath(data_sets.dirs[id], "ekipng")
    mkpath(root)

    # id, iN, iθ, icollumn)
    @infoe @sprintf("id = %d, iN = %d, iθ = %d, iz = %d", id1, 1, 1, ciz)
    z = get(data_sets, id1, 1, 1, ciz)

    local iI, iϵ, iκ, iIκ
    if idl == 0
        iI  = ciI
        iϵ  = ciϵ
        iκ  = ciκ
        iIκ = ciIκ
    elseif idl == 1
        iI  = ci1I
        iϵ  = ci1ϵ
        iκ  = ci1κ
        iIκ = ci1Iκ
    else
        iI  = ci2I
        iϵ  = ci2ϵ
        iκ  = ci2κ
        iIκ = ci2Iκ
    end

    iN, nθ = 1, 3
    #figure(1)
    II_10 = integrate_intensity(data_sets, id1, nθ, iN,   ciθ, iI, 400, "b", "r")
    II_11 = integrate_intensity(data_sets, id1, nθ, iN+1, ciθ, iI, 800, "g", "m")
    II_20 = integrate_intensity(data_sets, id2, nθ, iN,   ciθ, iI, 400, "b", "r")
    II_21 = integrate_intensity(data_sets, id2, nθ, iN+1, ciθ, iI, 800, "g", "m")
    #plot(II_10)
    #plot(II_11)
    #savefig(joinpath(root, "cubic_approximation.png"))
    #close(1)

    #out = open(joinpath(root, "result.txt"), "w")
    #write(out, @sprintf("I_400 - I_800 = %12.6e\n", II_0 - II_1))
    #@infoe @sprintf("with_emission = %s, I_400 - I_800 = %12.6e, %12.6e, %12.6e", with_emission, II_0 - II_1, II_0, II_1)

    for (j, zz) in enumerate(z)
        if zz > 0.2
            break
        end
    end

    j = 1
    zv = z[j:end]

    I11 = get(data_sets, id1, 1, 1, iI)[j:end]
    I12 = get(data_sets, id1, 2, 1, iI)[j:end]

    I21 = get(data_sets, id2, 1, 1, iI)[j:end]
    I22 = get(data_sets, id2, 2, 1, iI)[j:end]

    figure(2)
    loglog(z[j:end], I11,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], I12,  "g", label="800 ppm | ϵ>0")
    loglog(z[j:end], I21,  "m", label="400 ppm | ϵ=0")
    loglog(z[j:end], I22,  "y", label="800 ppm | ϵ=0")
    xlabel("h [m]")
    ylabel("I [W/m²]")
    legend(loc="best")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    savefig(joinpath(root, @sprintf("I_vs_h_%d.png", idl)))
    #0 : 10.0, 1: 6.4, 2: 10.0, 3: 10, 4: 10, 5: 6.37, 6: 10, 9: 9.12
    close(2)

    ϵ11  = get(data_sets, id1, 1, 1, iϵ)[j:end]
    ϵ12  = get(data_sets, id1, 2, 1, iϵ)[j:end]
    I1κ1 = get(data_sets, id1, 1, 1, iIκ)[j:end]
    I1κ2 = get(data_sets, id1, 2, 1, iIκ)[j:end]
    κ11  = get(data_sets, id1, 1, 1, iκ)[j:end]
    κ12  = get(data_sets, id1, 2, 1, iκ)[j:end]

    ϵ21  = get(data_sets, id2, 1, 1, iϵ)[j:end]
    ϵ22  = get(data_sets, id2, 2, 1, iϵ)[j:end]
    I2κ1 = get(data_sets, id2, 1, 1, iIκ)[j:end]
    I2κ2 = get(data_sets, id2, 2, 1, iIκ)[j:end]
    κ21  = get(data_sets, id2, 1, 1, iκ)[j:end]
    κ22  = get(data_sets, id2, 2, 1, iκ)[j:end]

    # ϵ Tv
    figure(3)
    semilogx(z[j:end], I1κ1-ϵ11,  "b", label="400 ppm | ϵ>0")
    #semilogx(z[j:end], I1κ2-ϵ12,  "g", label="800 ppm | ϵ>0")
    semilogx(z[j:end], I2κ1,  "m", label="400 ppm | ϵ=0")
    semilogx(z[j:end], I2κ2,  "y", label="800 ppm | ϵ=0")
    xlabel("h [m]")
    ylabel("(ϵ - κI) [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("eps-kapI_vs_h_%d.png",idl)))
    close(3)

    figure(4)
    plotfunc(z[j:end], ϵ11,  "b", label="400 ppm | ϵ>0")
    plotfunc(z[j:end], ϵ12,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("ϵ [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("eps_vs_h_%d.png", idl)))
    close(4)

    figure(5)
    loglog(z[j:end], κ11,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], κ12,  "g", label="800 ppm | ϵ>0")
    loglog(z[j:end], κ21,  "m", label="400 ppm | ϵ=0")
    loglog(z[j:end], κ22,  "y", label="800 ppm | ϵ=0")
    xlabel("h [m]")
    ylabel("κ [1/m]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("kap_vs_h_%d.png", idl)))
    close(5)

    figure(6)
    loglog(z[j:end], I1κ1,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], I1κ2,  "g", label="800 ppm | ϵ>0")
    loglog(z[j:end], I2κ1,  "m", label="400 ppm | ϵ=0")
    loglog(z[j:end], I2κ2,  "y", label="800 ppm | ϵ=0")
    xlabel("h [m]")
    ylabel("Iκ [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("Ikap_vs_h_%d.png", idl)))
    close(6)
end

function plot_results_single(data_sets, id, idl)
    with_emission = data_sets.parameter[id]["run_parameter"]["with_emission"]

    plotfunc = plot
    tickfunc = xticks

    root = joinpath(data_sets.dirs[id], "ekipng")
    mkpath(root)

    # id, iN, iθ, icollumn)
    @infoe @sprintf("id = %d, iN = %d, iθ = %d, iz = %d", id, 1, 1, ciz)
    z = get(data_sets, id, 1, 1, ciz)

    local iI, iϵ, iκ, iIκ
    if idl == 0
        iI  = ciI
        iϵ  = ciϵ
        iκ  = ciκ
        iIκ = ciIκ
    elseif idl == 1
        iI  = ci1I
        iϵ  = ci1ϵ
        iκ  = ci1κ
        iIκ = ci1Iκ
    else
        iI  = ci2I
        iϵ  = ci2ϵ
        iκ  = ci2κ
        iIκ = ci2Iκ
    end

    iN, nθ = 1, 3
    #figure(1)
    II_10 = integrate_intensity(data_sets, id, nθ, iN,   ciθ, iI, 400, "b", "r")
    II_11 = integrate_intensity(data_sets, id, nθ, iN+1, ciθ, iI, 800, "g", "m")
    #plot(II_10)
    #plot(II_11)
    #savefig(joinpath(root, "cubic_approximation.png"))
    #close(1)

    #out = open(joinpath(root, "result.txt"), "w")
    #write(out, @sprintf("I_400 - I_800 = %12.6e\n", II_0 - II_1))
    #@infoe @sprintf("with_emission = %s, I_400 - I_800 = %12.6e, %12.6e, %12.6e", with_emission, II_0 - II_1, II_0, II_1)

    for (j, zz) in enumerate(z)
        if zz > 0.2
            break
        end
    end

    j = 1
    zv = z[j:end]

    I11 = get(data_sets, id, 1, 1, iI)[j:end]
    I12 = get(data_sets, id, 2, 1, iI)[j:end]

    figure(2)
    loglog(z[j:end], I11,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], I12,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("I [W/m²]")
    legend(loc="best")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    savefig(joinpath(root, @sprintf("I_vs_h_%d.png", idl)))
    #0 : 10.0, 1: 6.4, 2: 10.0, 3: 10, 4: 10, 5: 6.37, 6: 10, 9: 9.12
    close(2)

    ϵ11  = get(data_sets, id, 1, 1, iϵ)[j:end]
    ϵ12  = get(data_sets, id, 2, 1, iϵ)[j:end]
    I1κ1 = get(data_sets, id, 1, 1, iIκ)[j:end]
    I1κ2 = get(data_sets, id, 2, 1, iIκ)[j:end]
    κ11  = get(data_sets, id, 1, 1, iκ)[j:end]
    κ12  = get(data_sets, id, 2, 1, iκ)[j:end]

    # ϵ Tv
    figure(3)
    semilogx(z[j:end], I1κ1-ϵ11,  "b", label="400 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("(ϵ - κI) [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("eps-kapI_vs_h_%d.png",idl)))
    close(3)

    figure(4)
    plotfunc(z[j:end], ϵ11,  "b", label="400 ppm | ϵ>0")
    plotfunc(z[j:end], ϵ12,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("ϵ [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("eps_vs_h_%d.png", idl)))
    close(4)

    figure(5)
    loglog(z[j:end], κ11,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], κ12,  "g", label="800 ppm | ϵ>0")

    xlabel("h [m]")
    ylabel("κ [1/m]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("kap_vs_h_%d.png", idl)))
    close(5)

    figure(6)
    loglog(z[j:end], I1κ1,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], I1κ2,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("Iκ [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("Ikap_vs_h_%d.png", idl)))
    close(6)
end

function run()
    #get_root_dirs()

    out_root = "/home/wester/ProjectsP/Julia/jllib/src/SimpleRadTrans/radoutput"
    subdir = "B"
    subsubdirs = readdir(joinpath(out_root, subdir))
    sort(subsubdirs)

    data_sets = makeDataSets(out_root, subdir, subsubdirs)

    k = collect(keys(data_sets.datasets))
    sort!(k)
    idl = 2
    for id in k
        plot_results_single(data_sets, 1, idl)
    end
end

run();
