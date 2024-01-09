using QuadGK
using JSON

import PyPlot as plt
function pyplot_gui(plt)
    plt.pygui(true)
    plt.pygui(:qt5)
end

using RWPhysConst

"""
    Read all result data from the log.log files
"""
mutable struct DataSets
    dirs      :: Vector{String}
    datasets  :: Dict
    parameter :: Dict
end
get_size(d::DataSets) = size(d,1)

"""
    load_results(datadir)

    read - result_outid_iN_iθ.hdf5
    read - input.json
"""
function load_results(datadir)
    files_and_dirs = readdir(datadir)
    data = []
    for f in files_and_dirs
        if occursin("result_", f)
            fs = split(splitext(f)[1], "_")
            iN = parse(Int64, fs[3])
            iθ = parse(Int64, fs[4])
            dat = load_array_as_hdf5(joinpath(datadir, f), group="A", dataset="A")
            #@infoe @sprintf("%d  %d  %s", iN, iθ, fs)
            push!(data, (iN, iθ, dat))
        end
    end

    parameter = JSON.parsefile(joinpath(datadir, "input.json"))#; dicttype=Dict, inttype=Int64, use_mmap=true)
    data, parameter
end

"""
    makeDataSets(out_root::String, subdir::String, subsubdirs::Vector{String})

    make datasets of all runs

    out_root   - root directory of data output
    subdir     - subdirectory
    subsubdirs - Vector of subsubdirectories

    returns - DataSets
"""
function makeDataSets(out_root::String, subdir::String, subsubdirs::Vector{String})
    dirs = []
    for ssd in subsubdirs
        push!(dirs, joinpath(out_root, subdir, ssd))
    end
    sort!(dirs)

    datasets = Dict()
    parameter = Dict()
    for (run_id, dd) in enumerate(dirs)
        datasets[run_id], parameter[run_id] = load_results(dd)
    end

    DataSets(dirs, datasets, parameter)
end

"""
    Get data column

    d - DataSets
    id {int} -- subsubdir
    iN {int} -- NCO2
    iθ {int} -- theta
    ie {int} -- column
        0     1  2  3  4  5  6  7  8    9
    iz, z, NCO2, θ, I, T, N, ϵ, κ, ΔλL, ΔλD

    array -- a column
"""
function get(d::DataSets, run_id, iN, iθ, ie)
    ds = d.datasets[run_id]

    #@infoe @sprintf("%d  %d  %d  %d  %d  %d  %d  %d", run_id, iN, iθ, ie, ds[1][1], ds[1][2], ds[2][1], ds[2][2])

    #ds = Vector{Tuple{Int64, Int64, Matrix}}
    i = -1
    for iv in eachindex(ds)
        #@infoe (ds[iv][1], ds[iv][2])
        if ds[iv][1] == iN && ds[iv][2] == iθ
            i = iv
        end
    end
    if i > -1
        # return : id, iN, iθ, ie =  z, NCO2, θ, I, T, N, ϵ, κ, ΔλL, ΔλD
        return ds[i][3][ie,:]
    end
    @warne ds[i]
end

# column indices
const c_z   =  1  # z
const c_C   =  2  # CO2 concentration
const c_θ   =  3  # θ
const c_T   =  4  # T
const c_N   =  5  # N
const c_I   =  8  # int_I
const c_ϵ   =  9  # int_ϵ
const c_κ   = 10  # int_κ
const c_Iκ  = 11  # int_Iκ
const c_I1  = 12  # int_I1
const c_ϵ1  = 13  # int_ϵ1
const c_κ1  = 14  # int_κ1
const c_Iκ1 = 15  # int_Iκ1
const c_I2  = 16  # int_I2
const c_ϵ2  = 17  # int_ϵ2
const c_κ2  = 18  # int_κ2
const c_Iκ2 = 19  # int_Iκ2
# integrated over all        wavelengths: int_I,  int_ϵ,  int_κ,  int_Iκ
# integrated over 1/6..n-1/6 wavelengths: int_I1, int_ϵ1, int_κ1, int_Iκ1
# integrated over 1/4..n-1/4 wavelengths: int_I2, int_ϵ2, int_κ2, int_Iκ2
"""
    Integrate intensity ove angle theta

    data_sets {[type]} -- [description]
    id        {[type]} -- [description]
    nθ        {[type]} -- [description]
    iN        {[type]} -- [description]
    NCO2      {[type]} -- [description]
    color1    {[type]} -- [description]
    color2    {[type]} -- [description]

    [type] -- [description]
"""
function integrate_intensity(data_sets::DataSets, run_id, iN, c_iI)
    θ_0 = get(data_sets, run_id, iN, 1, c_θ)[1]   # theta
    θ_1 = get(data_sets, run_id, iN, 2, c_θ)[1]   # theta
    θ_2 = get(data_sets, run_id, iN, 3, c_θ)[1]   # theta
    I_0 = get(data_sets, run_id, iN, 1, c_iI)[end] # intensity at TOA
    I_1 = get(data_sets, run_id, iN, 2, c_iI)[end] # intensity at TOA
    I_2 = get(data_sets, run_id, iN, 3, c_iI)[end] # intensity at TOA

    # qubic approximation of I(θ)
    R1 = I_1 - I_0
    R2 = I_2 - I_0

    a0  = I_0
    det = θ_1^2 * θ_2^3 - θ_1^3 * θ_2^2
    a2  = (R1 * θ_2^3 - R2 * θ_1^3) / det
    a3  = (R2 * θ_1^2 - R1 * θ_2^2) / det

    #x, w = gausslegendre(12)
    c1, err = quadgk(x -> cos(x)*sin(x)    , 0, π*0.5, rtol=1e-8)
    c2, err = quadgk(x -> cos(x)*sin(x)*x^2, 0, π*0.5, rtol=1e-8)
    c3, err = quadgk(x -> cos(x)*sin(x)*x^3, 0, π*0.5, rtol=1e-8)

    # integrated intensity
    Iint = 2.0*π * (a0*c1[1] + a2*c2[1] + a3*c3[1])

    return Iint
end

"""
    plot_results(data_sets::DataSets, run_id1, run_id2, id_λ_range) 

"""
function plot_results(data_sets::DataSets, run_id1, run_id2, id_λ_range)
    with_emission = data_sets.parameter[id]["run_parameter"]["with_emission"]

    plotfunc = plot
    tickfunc = xticks

    root = joinpath(data_sets.dirs[id], "ekipng")
    mkpath(root)

    # id, iN, iθ, icollumn)
    @infoe @sprintf("run_id = %d, iN = %d, iθ = %d, iz = %d", run_id1, 1, 1, ciz)
    z = get(data_sets, run_id1, 1, 1, ciz)

    local c_iI, c_iϵ, c_iκ, c_iIκ
    if id_λ_range == 0
        c_iI  = c_I
        c_iϵ  = c_ϵ
        c_iκ  = c_κ
        c_iIκ = c_Iκ
    elseif id_λ_range == 1
        c_iI  = c_I1
        c_iϵ  = c_ϵ1
        c_iκ  = c_κ1
        c_iIκ = c_Iκ1
    else
        c_iI  = c_I2
        c_iϵ  = c_ϵ2
        c_iκ  = c_κ2
        c_iIκ = c_Iκ2
    end

    iN = 1
    II_10 = integrate_intensity(data_sets, run_id1, iN,   c_iI)
    II_11 = integrate_intensity(data_sets, run_id1, iN+1, c_iI)

    II_20 = integrate_intensity(data_sets, run_id2, iN,   c_iI)
    II_21 = integrate_intensity(data_sets, run_id2, iN+1, c_iI)

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
    savefig(joinpath(root, @sprintf("I_vs_h_%d.png", id_λ_range)))
    #0 : 10.0, 1: 6.4, 2: 10.0, 3: 10, 4: 10, 5: 6.37, 6: 10, 9: 9.12
    close(2)

    #                              N, θ
    ϵ11  = get(data_sets, run_id1, 1, 1, c_iϵ)[j:end]
    ϵ12  = get(data_sets, run_id1, 2, 1, c_iϵ)[j:end]
    I1κ1 = get(data_sets, run_id1, 1, 1, c_iIκ)[j:end]
    I1κ2 = get(data_sets, run_id1, 2, 1, c_iIκ)[j:end]
    κ11  = get(data_sets, run_id1, 1, 1, c_iκ)[j:end]
    κ12  = get(data_sets, run_id1, 2, 1, c_iκ)[j:end]

    ϵ21  = get(data_sets, run_id2, 1, 1, c_iϵ)[j:end]
    ϵ22  = get(data_sets, run_id2, 2, 1, c_iϵ)[j:end]
    I2κ1 = get(data_sets, run_id2, 1, 1, c_iIκ)[j:end]
    I2κ2 = get(data_sets, run_id2, 2, 1, c_iIκ)[j:end]
    κ21  = get(data_sets, run_id2, 1, 1, c_iκ)[j:end]
    κ22  = get(data_sets, run_id2, 2, 1, c_iκ)[j:end]

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
    savefig(joinpath(root, @sprintf("eps-kapI_vs_h_%d.png",id_λ_range)))
    close(3)

    figure(4)
    plotfunc(z[j:end], ϵ11,  "b", label="400 ppm | ϵ>0")
    plotfunc(z[j:end], ϵ12,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("ϵ [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("eps_vs_h_%d.png", id_λ_range)))
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
    savefig(joinpath(root, @sprintf("kap_vs_h_%d.png", id_λ_range)))
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
    savefig(joinpath(root, @sprintf("Ikap_vs_h_%d.png", id_λ_range)))
    close(6)
end

function plot_results_single(data_sets, run_id, id_λ_range)
    with_emission = data_sets.parameter[run_id]["run_parameter"]["with_emission"]

    plotfunc = plot
    tickfunc = xticks

    root = joinpath(data_sets.dirs[run_id], "ekipng")
    mkpath(root)

    z = get(data_sets, run_id, 1, 1, c_z)

    local iI, iϵ, iκ, iIκ
    if id_λ_range == 0
        c_iI  = c_I
        c_iϵ  = c_ϵ
        c_iκ  = c_κ
        c_iIκ = c_Iκ
    elseif id_λ_range == 1
        c_iI  = c_I1
        c_iϵ  = c_ϵ1
        c_iκ  = c_κ1
        c_iIκ = c_Iκ1
    else
        c_iI  = c_I2
        c_iϵ  = c_ϵ2
        c_iκ  = c_κ2
        c_iIκ = c_Iκ2
    end

    iN = 1
    II_10 = integrate_intensity(data_sets, run_id, iN,   c_iI)
    II_11 = integrate_intensity(data_sets, run_id, iN+1, c_iI)

    for (j, zz) in enumerate(z)
        if zz > 0.2
            break
        end
    end

    j = 1
    zv = z[j:end]

    I11 = get(data_sets, run_id, 1, 1, c_iI)[j:end]
    I12 = get(data_sets, run_id, 2, 1, c_iI)[j:end]

    figure(2)
    loglog(z[j:end], I11,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], I12,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("I [W/m²]")
    legend(loc="best")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    savefig(joinpath(root, @sprintf("I_vs_h_%d.png", id_λ_range)))
    @infoe @sprintf("%e  %e  %e", I11[end], I12[end], I11[end]-I12[end])
    #0 : 10.0, 1: 6.4, 2: 10.0, 3: 10, 4: 10, 5: 6.37, 6: 10, 9: 9.12
    close(2)

    ϵ11  = get(data_sets, run_id, 1, 1, c_iϵ)[j:end]
    ϵ12  = get(data_sets, run_id, 2, 1, c_iϵ)[j:end]
    I1κ1 = get(data_sets, run_id, 1, 1, c_iIκ)[j:end]
    I1κ2 = get(data_sets, run_id, 2, 1, c_iIκ)[j:end]
    κ11  = get(data_sets, run_id, 1, 1, c_iκ)[j:end]
    κ12  = get(data_sets, run_id, 2, 1, c_iκ)[j:end]

    # ϵ Tv
    figure(3)
    semilogx(z[j:end], I1κ1-ϵ11,  "b", label="400 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("(ϵ - κI) [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("eps-kapI_vs_h_%d.png", id_λ_range)))
    close(3)

    figure(4)
    plotfunc(z[j:end], ϵ11,  "b", label="400 ppm | ϵ>0")
    plotfunc(z[j:end], ϵ12,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("ϵ [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("eps_vs_h_%d.png", id_λ_range)))
    close(4)

    figure(5)
    loglog(z[j:end], κ11,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], κ12,  "g", label="800 ppm | ϵ>0")

    xlabel("h [m]")
    ylabel("κ [1/m]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("kap_vs_h_%d.png", id_λ_range)))
    close(5)

    figure(6)
    loglog(z[j:end], I1κ1,  "b", label="400 ppm | ϵ>0")
    loglog(z[j:end], I1κ2,  "g", label="800 ppm | ϵ>0")
    xlabel("h [m]")
    ylabel("Iκ [W/m³]")
    #tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #title(dd[2])
    legend(loc="best")
    savefig(joinpath(root, @sprintf("Ikap_vs_h_%d.png", id_λ_range)))
    close(6)
end

function run(subdir)
    out_root = "/home/wester/Projects/GitHub/SARM/radoutput"
    subsubdirs = readdir(joinpath(out_root, subdir))
    sort(subsubdirs)

    data_sets = makeDataSets(out_root, subdir, subsubdirs)

    k = collect(keys(data_sets.datasets))
    sort!(k)
    id_λ_range = 0
    for run_id in k
        plot_results_single(data_sets, run_id, id_λ_range)
    end
end

subdir = "C"
run(subdir);
