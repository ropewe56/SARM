using PhysConst
using SimpleLog
using DataStructures

using PyPlot
pygui(true)
pygui(:qt5)

include("input/atmosphere.jl")
include("input/moleculardata.jl")
include("input/runparameter.jl")

include("profiles.jl")
include("resultdata.jl")
include("planck.jl")
include("spectrum.jl")

function run_radition_transfer()
    paths = OutPaths()
    
    rpar= RunParameter()
    
    atm = make_zpTN(rpar)

    mdH2O = get_TQ(rpar.H2Oiso, rpar.Tmin, rpar.Tmax, rpar.nT);
    mdCO2 = get_TQ(rpar.CO2iso, rpar.Tmin, rpar.Tmax, rpar.nT);
    
    H2O_line_data = get_line_data(rpar.H2Oout, rpar.λmin, rpar.λmax);
    CO2_line_data = get_line_data(rpar.CO2out, rpar.λmin, rpar.λmax);

    integrate(rpar, paths, atm, [mdH2O, mdCO2], [H2O_line_data, CO2_line_data])
end

