function make_outpaths(results)
    intensity = joinpath(results, "intensity")
    spectrum  = joinpath(results, "spectrum")
    atm       = joinpath(results, "atm")
    mkpath(results)
    mkpath(intensity)
    mkpath(spectrum)
    mkpath(atm)

    (
     results       =  results,
     intensity     = intensity,
     spectrum      = spectrum ,
     atm           = atm,
     logfile       = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/log.out",
     planck_single = joinpath(results,"planck_single.hdf5"),
     planck_multi  = joinpath(results,"planck_multi.hdf5")
    )
end
