include("preallocated.jl")
include("input/paths.jl")
include("input/parameter.jl")
include("input/atmosphere.jl")
include("input/moleculardata.jl")
include("input/mc_over_heigth.jl")
include("input/linedata.jl")

include("output/output.jl")

include("profiles.jl")
include("planck.jl")
include("spectrum.jl")
include("database.jl")

par = RunParameter()
par.θ         = deg2rad.([0.0]) # , 40.0, 80.0
par.planck_Ts = [288.0, 260.0, 240.0, 220.0]
par.species   = [:CO2]
par.c_ppm     = Dict(:H2O => fill(C0H2O_PPM, 3), )
par.c_ppm     = Dict(:CO2 => [430.0]) # 278.0, 430.0, 278.0*2.0

atm = Atmosphere(par);
datfiles = get_data_files()

mdH2O = MolecularData(:H2O, atm, datfiles[:H2O_Q], par.TQmin, par.TQmax);
mdCO2 = MolecularData(:CO2, atm, datfiles[:CO2_Q], par.TQmin, par.TQmax);

#hitran_to_hdf5(:H2O, datfiles[:H2Oout], datfiles[:H2Ohdf5], datfiles[:H2Ohdf5_compact], par.λmin, par.λmax, length(mdH2O.iso_a))
#hitran_to_hdf5(:CO2, datfiles[:CO2out], datfiles[:CO2hdf5], datfiles[:CO2hdf5_compact], par.λmin, par.λmax, length(mdCO2.iso_a))

H2O_line_data = LineData(datfiles[:H2Ohdf5_compact]);
CO2_line_data = LineData(datfiles[:CO2hdf5_compact]);

moleculardata = Dict(:H2O => mdH2O,         :CO2 => mdCO2)
linedata      = Dict(:H2O => H2O_line_data, :CO2 => CO2_line_data);

#moleculardata = Dict(:CO2 => mdCO2)
#linedata      = Dict(:CO2 => CO2_line_data);

outdir = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results"
set_paths!(par, outdir);
write_atm_to_hdf5(par.paths, atm)
rdb = create_results_db(par);
parameter_init(par)

integrate(par, rdb, atm, moleculardata, linedata)
