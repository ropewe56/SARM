include("include_sarm.jl")
include("postprocessing/plot_results.jl")

par = get_parameter()

par[:θ]         = deg2rad.([0.0]) # , 40.0, 80.0
par[:planck_Ts] = [288.0, 260.0, 240.0, 220.0, 215.0]
#par[:θ]         = deg2rad.(Vector{Float64}([0.0, 40.0, 80.0]))

par[:nh]        = 50

par[:concentrations] = Dict(:H2O => fill(C0H2O_PPM, 2), :CO2 => [400.0, 800.0])

par[:species]   = [:CO2]  #, :H2O]

par[:λmin]      = 12.0e-6
par[:λmax]      = 18.0e-6
par[:hmethod]   = :read
par[:hpath]     = "/home/wester/Projects/Julia/Climate-Energy/Sarm.rs/data/z.hdf5"
par[:h_iout]    = "/home/wester/Projects/Julia/Climate-Energy/Sarm.rs/data/z_iout.hdf5"

atmosphere      = Atmosphere(par);
molec_data_dict = get_molecular_data(par, atmosphere);

renew_hdf5 = true
line_data_dict  = get_line_data(par, molec_data_dict, renew_hdf5=renew_hdf5);

subdirs = readdir(OUTROOT)
#rm_subdirs(readdir(OUTROOT)[1:end])
#clear_subdir(subdir)
subdir = new_subdir_path()
subdir = subdirs[end]
par[:paths] = make_outpaths(subdir);

result_db = ResultDB(par);
parameter_init(par)

λb = make_λb(par)
create_planck_spectrum(par, λb);
Iλb = initial_intensity(par, λb);

write_atmosphere_to_hdf5(par[:paths], atmosphere);
prealloc = Preallocated();

integrate(par, result_db, λb, Iλb, atmosphere, molec_data_dict, line_data_dict, prealloc)

