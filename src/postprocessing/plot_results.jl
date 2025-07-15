import PyPlot as plt
plt.pygui(true)

pl = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results/2025-07-15T14:50/planck_single.hdf5"
plg = load_groups_as_hdf5(pl)
p = plg["intensity"]
pλ, pI = p["wl"], p["I"]
plt.plot(pλ, pI)

dbpath = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results/2025-07-15T14:50/db.sqlite3"
db = open_db(dbpath)
for n in names(df) println(n) end

df = select_from_rdb(db, ic=2, iθ=1)
h = df[!,"h"]
IH2O = df[!,"int_IH2O"]
ICO2 = df[!,"int_ICO2"]

plt.plot(h,I)
plt.plot(h,df[!,"T"])

