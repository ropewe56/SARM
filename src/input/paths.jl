using Dates
using Printf

function make_outpaths(outdir)
    d = Dates.now()
    subdir = @sprintf("%s", d)[1:16]

    root      = joinpath(outdir, subdir)
    intensity = joinpath(root, "intensity")
    spectrum  = joinpath(root, "spectrum")
    atm       = joinpath(root, "atm")
    mkpath(root)
    mkpath(intensity)
    mkpath(spectrum)
    mkpath(atm)

    (outroot           = root,
     intensity         = intensity,
     spectrum          = spectrum,
     atm               = atm,
     logfile           = joinpath(root, "log.out"),
     dbpath            = joinpath(root, "db.sqlite3"),
     planck_single     = joinpath(intensity, "planck_single.hdf5"),
     planck_multi      = joinpath(intensity, "planck_multi.hdf5"),
     initial_intensity = joinpath(intensity, "initial_intensity.hdf5")
    )
end
