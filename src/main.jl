include("include_sarm.jl")
include("postprocessing/plot_results.jl")

par = get_parameter()

par[:θ]         = deg2rad.([0.0]) # , 40.0, 80.0
par[:planck_Ts] = [288.0, 260.0, 240.0, 220.0, 215.0]
#par[:θ]         = deg2rad.(Vector{Float64}([0.0, 40.0, 80.0]))

par[:nh]        = 50

par[:concentrations] = Dict(:H2O => fill(C0H2O_PPM, 3), :CO2 => [278.0, 430.0, 278.0*2.0])

par[:species]   = [:CO2]  #, :H2O]

par[:λmin]      = 13.5e-6
par[:λmax]      = 16.5e-6

atmosphere      = Atmosphere(par);
molec_data_dict = get_molecular_data(par, atmosphere);

#line_data = get_species_line_data(par, :CO2, molec_data_dict[:CO2]; renew_hdf5=true);
#line_data = get_species_line_data(par, :H2O, molec_data_dict[:H2O]; renew_hdf5=true);
line_data_dict  = get_line_data(par, molec_data_dict);

subdirs = readdir(OUTROOT)
#rm_subdirs(readdir(OUTROOT)[1:end])
#clear_subdir(subdir)
subdir = new_subdir_path()
par[:paths] = make_outpaths(subdir);

result_db = ResultDB(par);
parameter_init(par)

result_db.colnames

λb = make_λb(par)
create_planck_spectrum(par, λb);
Iλb = initial_intensity(par, λb);

write_atmosphere_to_hdf5(par[:paths], atmosphere);
prealloc = Preallocated();

integrate(par, result_db, λb, Iλb, atmosphere, molec_data_dict, line_data_dict, prealloc)

