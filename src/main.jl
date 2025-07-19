include("include_sarm.jl")

par = get_parameter()

par[:θ]         = deg2rad.([0.0]) # , 40.0, 80.0
par[:planck_Ts] = [288.0, 260.0, 240.0, 220.0]
#par[:θ]         = deg2rad.(Vector{Float64}([0.0, 40.0, 80.0]))

par[:species]   = [:CO2]
par[:c_ppm]     = Dict(:H2O => fill(C0H2O_PPM, 3), )
par[:c_ppm]     = Dict(:CO2 => [430.0]) # 278.0, 430.0, 278.0*2.0

atm             = Atmosphere(par);
molec_data_dict = get_molecular_data(par);
line_data_dict  = get_line_data(par, molecular_data);

par[:paths] = make_outpaths();
write_atm_to_hdf5(par[:paths], atm)
rdb = create_results_db(par);
parameter_init(par)


integrate(par, rdb, atm, molec_data_dict, line_data_dict)
