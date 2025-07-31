using SimpleLog
using SpecialFileIO
using Printf
using DataFrames

import PyPlot as plt
plt.pygui(true)

function get_list()
    readdir(OUTROOT)
end

root = get_list()
iroot = joinpath(root[end], "spectrum")

dbpath = joinpath(OUTROOT, root[end], "db.sqlite3")
db = open_db(dbpath)

function sort_hdf5_paths(h5p)
    nb = []
    for f in h5p1
        push!(nb, parse(Int, split(basename(f), '_')[4]))
    end
    index = sortperm(nb)
    h5p[index]
end

function get_dbresults(db, ic)
    df = select_from_rdb(db, ic=ic, iθ=1)
    h      = df[!,"h"]
    int_I  = df[!,"int_I"]
    int_ϵ  = df[!,"int_ϵ"]
    int_Iκ = df[!,"int_Iκ"]
    h5p    = filter(x -> x != "none", df[!,"hdf5_path"])
    h5p, h, int_I, int_ϵ, int_Iκ
end

h5p1, h1, I1, ϵ1, κ1 = get_dbresults(db, 1)
h5p2, h2, I2, ϵ2, κ2 = get_dbresults(db, 2)

h5p1 = sort_hdf5_paths(h5p1)
h5p2 = sort_hdf5_paths(h5p2)

plt.plot(h1, I1)
plt.plot(h2, I2)

####
f0 = load_groups_as_hdf5(h5p1[end])
λ0 = f0["sarm"]["λ"]
I0 = f0["sarm"]["I"]

f1 = load_groups_as_hdf5(h5p2[end])
λ1 = f1["sarm"]["λ"]
I1 = f1["sarm"]["I"]

plt.plot(λ0, I0)
plt.plot(λ1, I1)
plt.show()

