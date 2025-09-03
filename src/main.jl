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

    projroot = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl"
    dbpath = joinpath(projroot, "data", "MolecularData.db")
    db = SQLite.DB(dbpath)

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
    molec_data_dict = get_molecular_data(db, par, atmosphere);
    # save input data
    save_input_to_hdf5(par.p.input_data, atmosphere, molec_data_dict)

    # line data
    line_data_dict = get_dbline_data(db, par.m.species, par.w.λmin, par.w.λmax);

    result_db = ResultDB(par);
    parameter_init_and_save(par);

    λb = make_λb(par)
    create_planck_spectrum(par, λb);
    Iλb0 = initial_intensity(par, λb);

    subdir, par, result_db, λb, Iλb0, atmosphere, molec_data_dict, line_data_dict
end

species = [:CO2,:H2O]
jl_rs = "jl"

subdir, par, result_db, λb, Iλb0, atmosphere, molec_data_dict, line_data_dict = init_sarm(species, jl_rs);

integrate(par, result_db, λb, Iλb0, atmosphere, molec_data_dict, line_data_dict)

#show_all_par(par)
#start_rust(subdir)

#julia --optimize=3 --inline=yes --check-bounds=no --math-mode=fast --threads=10 src/main.jl

E2 = 1.3218888896425107e-20

E1 = 7.3708e-20
E2 = 8.6920e-20

exp(-E1/(c_kB*287))
exp(-E2/(c_kB*287))

@infoe iso, λ21, A21, ϵ, κ, E1, E2, Niso, N, N1, N2, Qiso[iso], miso[iso], aiso[iso], cspech
@infoe T, p, exp(- E1 * β), exp(- E2 * β)

(1, 1.5035712033956973e-5, 1.591, 0.054617477364122737, 9.420461204727894e-9, 7.37084600687193e-24, 1.3218888896425107e-20, 
1.0028868725110011e22, 2.5474992673541107e25, 8.33112194613359e20, 3.265266170130596e19, 
276.3574824716923, 7.304683127333609e-26, 0.984204, 0.0003999933333333333)
(287.9972230769231, 101294.45000000001, 0.9981479901783495, 0.03599132716207241)

1,  1.50357124e-5,   1.5910e0, 4.8697e-10, 8.3993e-17, 7.3708e-20, 8.6920e-20,  1.0029e22,  7.4281e12,  2.9113e11,   2.7636e2, 7.3047e-26,  9.8420e-1,  3.9999e-4
 2.88e2,   1.01e5,  8.90e-9, 3.21e-10

