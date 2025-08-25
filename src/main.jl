include("include_sarm.jl")
include("postprocessing/plot_results.jl")
include("start_rs.jl")

function subdir_util()
    subdirs = readdir(OUTROOT)
    rm_subdirs(readdir(OUTROOT)[1:end])
    clear_subdir(subdirs[end])
    subdir = subdirs[end]
    subdir
end

function init_sarm(species, jl_rs)
    result_root = if jl_rs == "jl"
        RESULT_JL
    else
        RESULT_RS
    end

    par = SarmParameter();

    par.r.θ         = deg2rad.([0.0]) # , 40.0, 80.0
    #par[:θ]         = deg2rad.(Vector{Float64}([0.0, 40.0, 80.0]))
    par.c.planck_Ts = [288.0, 260.0, 240.0, 220.0, 215.0]

    par.h.nh     = 500
    par.h.hmethod = :dh

    par.w.λmin    = 12.0e-6
    par.w.λmax    = 18.0e-6
    par.m.c_ppm   = Dict(:H2O => fill(C0H2O_PPM, 2).*PPM, :CO2 => [400.0, 800.0].*PPM)
    par.m.species = species

    subdir = @sprintf("%s_%s", jl_rs, join(par.m.species, "_"))
    par.p = make_outpaths(result_root, subdir);

    # input data
    atmosphere      = Atmosphere(par.h);
    molec_data_dict = get_molecular_data(par.m, atmosphere);
    # save input data
    save_input_to_hdf5(par.p.input_data, atmosphere, molec_data_dict)

    # line data
    renew_hdf5 = false#true
    line_data_dict  = get_line_data(par, molec_data_dict, renew_hdf5=renew_hdf5);

    result_db = ResultDB(par);
    parameter_init_and_save(par);

    λb = make_λb(par)
    create_planck_spectrum(par, λb);
    Iλb0 = initial_intensity(par, λb);

    subdir, par, result_db, λb, Iλb0, atmosphere, molec_data_dict, line_data_dict
end

species = [:CO2,:H2O]
jl_rs = "rs"

subdir, par, result_db, λb, Iλb0, atmosphere, molec_data_dict, line_data_dict = init_sarm(species, jl_rs);

integrate(par, result_db, λb, Iλb0, atmosphere, molec_data_dict, line_data_dict)

#show_all_par(par)
#start_rust(subdir)

