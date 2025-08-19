using SimpleLog
using SpecialFileIO
using Printf
using DataFrames

import PyPlot as plt
plt.pygui(true)

include("../include_sarm.jl")

function get_result_root(jl_rs, species)
    result_root = if jl_rs == "jl"
        RESULT_JL
    else
        RESULT_RS
    end
    dirs = readdir(result_root)

    subdir = @sprintf("%s_%s", jl_rs, join(species, "_"))
    result_subdir = nothing
    for d in dirs
        if d == subdir
            result_subdir = joinpath(result_root, d)
        end
    end
    result_subdir
end

function sort_hdf5_paths(h5p)
    nb = []
    for f in h5p
        push!(nb, parse(Int, split(basename(f), '_')[4]))
    end
    index = sortperm(nb)
    h5p[index]
end

function load_database(result_root)
    dbpath = joinpath(result_root, subdir, "db.sqlite3")
    db = open_db(dbpath)
    db
end

function get_dbresults(db; ic=1, iθ=1)
    df = select_from_rdb(db; ic=ic, iθ=iθ)
    h      = df[!,"h"]
    int_I  = df[!,"int_I"]
    int_ϵ  = df[!,"int_ϵ"]
    int_Iκ = df[!,"int_Iκ"]
    h5p    = filter(x -> x != "none", df[!,"hdf5_path"])
    h5p, h, int_I, int_ϵ, int_Iκ
end

function plot_I_vs_h(db, ic_label)
    results = []
    for (ic, iθ, lab) in ic_label
        h5p1, h1, I1, ϵ1, κ1 = get_dbresults(db, ic=ic, iθ=iθ);
        plt.plot(h1, I1, label=lab)
        push!(results, (h5p1, h1, I1, ϵ1, κ1))
    end
    plt.legend()

    results
end

species = [:CO2,:H2O]

jl_rs = "jl"
result_root1 = get_result_root(jl_rs, species)
dbpath1 = joinpath(result_root1, "db.sqlite3")
db1 = open_db(dbpath1)
db = db1
ic=1; iθ = 1; label="jl CO2 = 400 ppm"

results1 = plot_I_vs_h(db1, ( (ic=1, iθ = 1, label="jl CO2 = 400 ppm"), 
                              (ic=2, iθ = 1, label="jl CO2 = 800 ppm")));
results1[1][3][end]
results1[2][3][end]


jl_rs = "rs"
result_root2 = get_result_root(jl_rs, species)
dbpath2 = joinpath(result_root2, "db.sqlite3")
db2 = open_db(dbpath2)
results2 = plot_I_vs_h(db2, ( (ic=1, iθ = 1, label="rs CO2 = 400 ppm"), 
                              (ic=2, iθ = 1, label="rs CO2 = 800 ppm")));
results2[1][3][end]
results2[2][3][end]

plt.xlabel("h [m]")
plt.ylabel("I [W/m^2]")

plt.savefig("I_vs_h.png")

#CO2
#27.726757150765113
#26.649775087752026

#CO2,H2O
#27.60868716167992
#26.564405756108336
##===================

pl = load_groups_as_hdf5("/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results/jl_CO2_H2O/intensity/planck_multi.hdf5")

Tp = pl["TλI"]["T"]
λp = pl["TλI"]["λ"]
Ip = pl["TλI"]["I"]
n1, n2 = size(Ip)
i = 2
for i in 1:n2
    Ipi = Ip[:,i]
    Ti = @sprintf("T = %6.1f", Tp[i])
    plt.plot(λp, Ipi, label=Ti)
end

# h5p, h, int_I, int_ϵ, int_Iκ
h5p1 = results1[1][1]
h5p2 = results2[1][1]
h5p1 = sort_hdf5_paths(h5p1);
h5p2 = sort_hdf5_paths(h5p2);

f1 = load_groups_as_hdf5(h5p1[end]);
f2 = load_groups_as_hdf5(h5p2[end]);

#keys(f1["sarm"])
#show(IOContext(stdout, :limit=>false), MIME"text/plain"(), keys(f1["sarm"]))
#show(IOContext(stdout, :limit=>false), MIME"text/plain"(), keys(f2["sarm"]))


λb1 = f1["sarm"]["λb"]
λb2 = f2["sarm"]["λb"]
Ib1 = f1["sarm"]["Iλb"]
Ib2 = f2["sarm"]["Iλb"]

plt.plot(λb1, Ib1)
plt.plot(λb2, Ib2)
plt.axis([nothing, nothing, 0.0, nothing])

extrema(λl1 .- λl2)

ϵl1 = f1["sarm"]["ϵl_CO2"]
ϵl2 = f2["sarm"]["ϵl_CO2"]
extrema(ϵl1.-ϵl2)

κl1 = f1["sarm"]["κl_CO2"]
κl2 = f2["sarm"]["κl_CO2"]
extrema(κl1.-κl2)

infl1 = f1["sarm"]["intf_CO2"]
infl2 = f2["sarm"]["intf_CO2"]
length(infl1)
length(infl2)
extrema(infl1.-infl2)

λb1 = f1["sarm"]["λb"]
ϵb1 = f1["sarm"]["ϵb_CO2"]
κb1 = f1["sarm"]["κb_CO2"]

λb2 = f2["sarm"]["λb"]
ϵb2 = f2["sarm"]["ϵb_CO2"]
κb2 = f2["sarm"]["κb_CO2"]

extrema(λb1.-λb2)
extrema(ϵb1.-ϵb2)
extrema(κb1.-κb2)

i = argmax(ϵb1.-ϵb2)
ϵb1[i]
ϵb2[i]

plt.plot(λb1, ϵb1)
plt.plot(λb2, ϵb2)

