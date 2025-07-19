using SimpleLog
using SpecialFileIO
using Printf
import PyPlot as plt
plt.pygui(true)

results_root = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results"

function get_list()
    readdir(results_root)
end

function plot_planck(root, λ1, λ2)
    path_pls = joinpath(results_root, root, "intensity", "planck_single.hdf5")
    path_plm = joinpath(results_root, root, "intensity", "planck_multi.hdf5")
    path_Ii  = joinpath(results_root, root, "intensity", "initial_intensity.hdf5")

    pls = load_groups_as_hdf5(path_pls)
    plm = load_groups_as_hdf5(path_plm)
    pli = load_groups_as_hdf5(path_Ii)

    p = pls["TλI"]
    λs, Is, Ts = p["λ"], p["I"], p["T"]

    p = plm["TλI"]
    λm, Im, Tm = p["λ"], p["I"], p["T"]

    p = pli["TλI"]
    λi, Ii, Ti = p["λ"], p["I"], p["T"]

    plt.figure()
    plt.plot(λs, Is, label = "Is")
    for (i,T) in enumerate(Tm)
        plt.plot(λm, Im[:,i], label = @sprintf("T = %5.1f", Tm[i]))
    end
    plt.plot(λi, Ii, label = "Ii")
    if λ1 > 0.0 && λ2 > 0.0
        plt.xlim(λ1, λ2)
    end
    plt.legend()
end

function get_results(root, ic, iθ, hi)
    dbpath = joinpath(results_root, root, "db.sqlite3")
    db = open_db(dbpath)
    df = select_from_rdb(db, ic=1, iθ=1)

    for n in names(df)
        @printf("%s  ", n)
    end
    @printf("\n")

    h  = df[!,"h"]
    ih = max(1, min(length(h), argmin(abs.(h .- hi))))

    species = split(df[1,"species"], ",")
    for spec in species
        cih = df[!,"cih"*spec][ih]
        I   = df[!,"int_I"*spec][ih]
        κ   = df[!,"int_κ"*spec][ih]
        ϵ   = df[!,"int_ϵ"*spec][ih]
        Iκ  = df[!,"int_Iκ"*spec][ih]
        ΔλD = df[!,"ΔλD"*spec][ih]
        ΔλL = df[!,"ΔλL"*spec][ih]
        @infoe spec, cih, I, κ, ϵ, Iκ, ΔλL, ΔλD
    end

    species, df[!,"hdf5_path"], h, ih
end

function plot_result(hdf5_path)
    path_Ii  = joinpath(results_root, root, "intensity", "initial_intensity.hdf5")
    pli = load_groups_as_hdf5(path_Ii)
    p = pli["TλI"]
    λi, Ii = p["λ"], p["I"]

    groups = load_groups_as_hdf5(hdf5_path)
    data = groups["sarm"]
    keys(data)

    λ = data["λ"]
    I = data["I"]
    ϵ = data["ϵ"]
    κ = data["κ"]
    ϵ_CO2 = data["ϵ_CO2"]
    κ_CO2 = data["κ_CO2"]

    λl = data["λl_CO2"]
    Sl = data["Sl_CO2"]
    ϵl = data["ϵl_CO2"]
    κl1 = data["κl1_CO2"]
    κl2 = data["κl2_CO2"]

    plt.figure()

    plt.plot(λ,I)
    plt.plot(λi,Ii)
    plt.xlabel("λ")
    plt.ylabel("I")

    plt.figure()
    plt.plot(λ,ϵ)
    plt.xlabel("λ")
    plt.ylabel("ϵ")

    plt.figure()
    plt.plot(λ,κ.*I)
    plt.xlabel("λ")
    plt.ylabel("κI")

    plt.figure()
    plt.plot(λl,ϵl)
    plt.xlabel("λ")
    plt.ylabel("ϵl")

    plt.figure()
    plt.plot(λl,κl1)
    plt.xlabel("λ")
    plt.ylabel("κl1")

    plt.figure()
    plt.plot(λl,κl2)
    plt.xlabel("λ")
    plt.ylabel("κl2")

    plt.figure()
    plt.plot(λl,Sl)
    plt.xlabel("λ")
    plt.ylabel("Sl")
end

root = readdir(results_root)[end]

#plot_planck(root, 10.0e-6, 20.0e-6)

species, hdf5_paths, h, ih = get_results(root, 1, 1, 70000.0);
hdf5_path = hdf5_paths[ih]
plot_result(hdf5_path)
