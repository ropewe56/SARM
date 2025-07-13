using Parameters

@with_kw mutable struct OutPaths
    datadir = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/data"

    H2Oiso = joinpath(datadir, "H2O", "H2O_Q", "H2O_Isotopes.txt")
    CO2iso = joinpath(datadir, "CO2", "CO2_Q", "CO2_Isotopes.txt")
    H2Oout = joinpath(datadir, "H2O", "H2O_rwfmt_ISO-0-1.out")
    CO2out = joinpath(datadir, "CO2", "CO2_rwfmt_ISO-0-12_wl-12-18-mum.out")

    results   = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results"
    intensity = joinpath(results, "intensity")
    spectrum  = joinpath(results, "spectrum")
    atm       = joinpath(results, "atm")

    logfile       = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/log.out"
    planck_single = joinpath(results,"planck_single.hdf5")
    planck_multi  = joinpath(results,"planck_multi.hdf5")
end

function OutPaths(results)
    intensity = joinpath(results, "intensity")
    spectrum  = joinpath(results, "spectrum")
    atm       = joinpath(results, "atm")
    paths     = OutPaths()
    paths.results   = results
    paths.intensity = intensity
    paths.spectrum  = spectrum 
    paths.atm       = atm      
    mkpath(paths.results)
    mkpath(paths.intensity)
    mkpath(paths.spectrum)
    mkpath(paths.atm)
    paths
end
