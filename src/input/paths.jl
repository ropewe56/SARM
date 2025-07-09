using Parameters

@with_kw mutable struct OutPaths
    datadir = "/home/wester/Projects/Julia/Climate-Energy/SARM/data"
    H2Oiso = joinpath(datadir, "H2O", "H2O_Q", "H2O_Isotopes.txt")
    CO2iso = joinpath(datadir, "CO2", "CO2_Q", "CO2_Isotopes.txt")
    H2Oout = joinpath(datadir, "H2O", "H2O_rwfmt_ISO-0-1.out")
    CO2out = joinpath(datadir, "H2O", "CO2_rwfmt_ISO-0-12_wl-12-18-mum.out")

    outdir = "/home/wester/Projects/Julia/Climate-Energy/SARM/radoutput"
    intensity_dir = joinpath(outdir, "intensity")
    planck_single = "planck_single.hdf5"
    planck_multi = "planck_multi.hdf5"
end