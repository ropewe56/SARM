include("preallocated.jl")
include("input/paths.jl")
include("input/parameter.jl")
include("input/atmosphere.jl")
include("input/moleculardata.jl")
include("input/mc_over_heigth.jl")
include("input/linedata.jl")

include("output/output.jl")
include("resultdata.jl")

include("profiles.jl")

include("planck.jl")
include("spectrum.jl")

outdir = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results1"
par = make_parameter(outdir);

par.θ         = deg2rad.([0.0, 40.0, 80.0])
par.planck_Ts = [288.0, 260.0, 240.0, 220.0]
par.species   = [:H2O, :CO2]
par.c_ppm     = Dict(:H2O => fill(C0H2O_PPM, 3), :CO2 => [278.0, 430.0, 278.0*2.0])
par.c_ppm     = Dict(:CO2 => [278.0, 430.0, 278.0*2.0])

atm = Atmosphere(par);
write_atm_to_hdf5(par.paths, atm)

mdH2O = MolecularData(:H2O, atm, H2O_Q, par.TQmin, par.TQmax);
mdCO2 = MolecularData(:CO2, atm, CO2_Q, par.TQmin, par.TQmax);

H2O_line_data = LineData(:H2O, H2Oout, par.λmin, par.λmax, length(mdH2O.iso_a));
CO2_line_data = LineData(:CO2, CO2out, par.λmin, par.λmax, length(mdCO2.iso_a));

moleculardata = Dict(:H2O => mdH2O,         :CO2 => mdCO2)
linedata      = Dict(:H2O => H2O_line_data, :CO2 => CO2_line_data);

moleculardata = Dict(:CO2 => mdCO2)
linedata      = Dict(:CO2 => CO2_line_data);

integrate(par, atm, moleculardata, linedata)


hdf5_path = SpecialFileIO.make_hdf5_path(path)
path = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results1/spectrum/spectrum_001_001_1_1_ 0.0.hdf5"
fr = load_groups_as_hdf5(path)
keys(fr["I"])

import PyPlot as plt
plt.pygui(true)

wl = fr["I"]["wl"]

I = fr["I"]["e_CO2"]
I = fr["I"]["I"]

plt.plot(wl, I)