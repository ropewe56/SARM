using SimpleLog

using PyPlot
pygui(true)
pygui(:qt5)

include("input/paths.jl")
include("input/runparameter.jl")
include("input/atmosphere.jl")
include("input/moleculardata.jl")
include("input/linedata.jl")

include("profiles.jl")
include("resultdata.jl")

include("planck.jl")
include("spectrum.jl")


function run_radition_transfer()
    paths = OutPaths()
    par= RunParameter()
    
    atm = Atmosphere(par)

    mdH2O = MolecularData(paths.H2Oiso, par.TQmin, par.TQmax)
    mdCO2 = MolecularData(paths.CO2iso, par.TQmin, par.TQmax);
    
    H2O_line_data = LineData(paths.H2Oout, par.λmin, par.λmax);
    CO2_line_data = LineData(paths.CO2out, par.λmin, par.λmax);

    md = [mdH2O, mdCO2];
    ld = [H2O_line_data, CO2_line_data];
    
    integrate(par, paths, atm, md, ld)
end
