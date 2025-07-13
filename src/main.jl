include("preallocated.jl")
include("input/paths.jl")
include("input/parameter.jl")
include("output/output.jl")
include("input/atmosphere.jl")
include("input/moleculardata.jl")
include("input/linedata.jl")

include("profiles.jl")
include("resultdata.jl")

include("planck.jl")
include("spectrum.jl")

paths = OutPaths("/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results1");

par = RunParameter();
par.θ = deg2rad.([0.0, 40.0, 80.0])
par.planck_Ts = [288.0, 260.0, 240.0, 220.0]

par.c_ppm = set_concentrations([[ch0H2O_ppm], [278.0, 430.0, 278.0*2.0]])
cCO2_ppm = [278.0, 430.0, 278.0*2.0]
cH2O_ppm = fill(ch0H2O_ppm, 3)

atm = Atmosphere(par);

mdH2O = MolecularData(:H2O, cH2O_ppm, atm, paths.H2Oiso, par.TQmin, par.TQmax);
mdCO2 = MolecularData(:CO2, cCO2_ppm, atm, paths.CO2iso, par.TQmin, par.TQmax);

#par.c_ppm = set_concentrations([[278.0, 430.0, 278.0*2.0]])

atm = Atmosphere(par);

mdH2O = MolecularData(:H2O, atm, paths.H2Oiso, par.TQmin, par.TQmax, atm.hi);
mdCO2 = MolecularData(:CO2, atm, paths.CO2iso, par.TQmin, par.TQmax, atm.hi);

H2O_line_data = LineData(paths.H2Oout, par.λmin, par.λmax, length(mdH2O.iso_a));
CO2_line_data = LineData(paths.CO2out, par.λmin, par.λmax, length(mdCO2.iso_a));

moleculardata = [mdH2O, mdCO2];
linedata = [H2O_line_data, CO2_line_data];

#moleculardata = [mdCO2];
#linedata = [CO2_line_data];

integrate(par, paths, atm, moleculardata, linedata)


