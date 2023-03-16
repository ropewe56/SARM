using Common
using PhysConst
using PyPlot
pygui(true)
pygui(:qt5)

include("create_input/create_QTp_lines_files.jl")
include("create_input/create_Nthetaz_json_files.jl")

include("interpolator.jl")
include("parameter.jl")
include("profiles.jl")
include("resultdata.jl")
include("spectralline.jl")
include("planck.jl")
include("spectrum.jl")

#path1 = "/home/wester/Projects/Julia/Private/SimpleRadTrans.jl/HITRAN/CO2_rwfmt_ISO-0-12_wl-12-18-mum.npy"
#isfile(path1)
#data = read_npz(path1; set_extension=false)
#
#path2 = "/home/wester/Projects/Julia/Private/SimpleRadTrans.jl/HITRAN/CO2_rwfmt_ISO-0-12_wl-12-18-mum.hdf5"
#write_hdf5(path2, Dict("HITRAN" => Dict("CO2" => data)))
#
#grdst = read_hdf5(path2)
#print_names_and_types(grdst)


"""
    plotz()
plot p and T as fucntion of height z
"""
function plotz()
    data_dir, input_dir, out_root = get_paths()
    hp = read_npy(joinpath(input_dir, "h_p.npy"))
    hT = read_npy(joinpath(input_dir, "h_T.npy"))
    println(size(hp))

    zmin, zmax, dzmin, dzmax, n, e = 0.0, 7.0e4, 10.0, 1.0e3, 200, 3.0

    h = collect(LinRange(zmin, zmax, n))

    plot(hp[1,:], hp[2,:])
    figure()
    plot(hT[1,:], hT[2,:])
    figure()

    z, z_iout, zh = make_z(zmin, zmax, dzmin, dzmax, n, e)

    plot(h, zh, "b", label="exp")

    plot(h, z, "r", label="z")
    legend()
end
#plotz()

function make_histograms()
    data_dir, input_dir, out_root = get_paths()

    # parameter
    max_isotope_id = 11
    ln_bins        = 80
    En_bins        = 100

    # read line data
    specdat_path = joinpath(data_dir, "CO2_rwfmt_ISO-0-12_wl-12-18-mum.npy")
    ld = LineData(specdat_path, max_isotope_id)

    # write
    hdf5_path = joinpath(input_dir, "EE_histogram.hdf5")
    edges, weights = energy_energy_histogram(ld, ln_bins, En_bins, hdf5_path)

    fig, ax1 = plt.subplots(figsize=(8, 6), ncols=1)
    pos = ax1.imshow(weights, origin="lower", interpolation ="none", aspect="auto",  extent=[12.0, 18.0, 0.0, edges[2][end]/c_e])
    fig.colorbar(pos, ax=ax1)
    ax1.set_xlabel("transition wavelength λ [μm]")
    ax1.set_ylabel("energy of lower level [eV]")
    ax1.set_title("number of lines")
    fig.savefig(joinpath(input_dir, "EE_histogram.png"))
end
#make_histograms()

"""
Input parameter:
    CO2 concentration values
    Inclination angle values
    Input parameter for determining height values
        T_of_h = false: temperature and density values vary according to height in the atmosphere
               = true : temperature is constant, density value varies according to height in the atmosphere
        with_emission = true: include emission in the radiation transfer equation
                        false: neglect emission in the radiation transfer equation
        albedo:
        nisos: number of CO2 isotopes included
        bg: background added to the absorption and emission coefficients because of the cut off of Lorentz shape contributions
        Δλ: wavelength resolution
"""
input_run_parameter = Dict(
    "fixed_parameter" => Dict(
        "λmin"      => 1.2e-5,
        "λmax"      => 1.8e-5,
        "Δλ_factor" => 20.0,       # relative width of the line shapes
        "ΔλL"       => 1.0e-9,     # Lorentz width
        "p_ref"     => 1013.25e2,  # N / m^2
        "T_ref"     => 296.0,      # K
        "CO2_mass"  => (12.011 + 2.0 * 15.999) * 1.660539040e-27,
        "integrate" => true,
        "T_surface" => 288.0,
        "Planck_Ts" => [200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 250.0, 260.0, 270.0, 280.0, 288.0],
    ),
    "npy_parameter" => Dict(
        "NCO2"      => [400.0*1.0e-6, 800.0*1.0e-6],
        "theta_deg" => [0.0, 40.0, 80.0],
        "z"         => Dict("zmin"=> 0.0, "zmax"=> 7.0e4, "dzmin"=> 10.0, "dzmax"=> 1.0e3, "n" => 200,  "exponent" => 1),
    ),
    "loop_parameter" => [
        Dict("T_of_h" => true,  "N_of_h" => true, "with_emission" => true,  "albedo" => 0.0, "max_isotope_id" =>  11, "back_ground"=> 0.0, "Δλ" => 1.0e-11, "initial_intensity" => "zero"),
        Dict("T_of_h" => true,  "N_of_h" => true, "with_emission" => true,  "albedo" => 0.0, "max_isotope_id" =>  11, "back_ground"=> 0.0, "Δλ" => 1.0e-11, "initial_intensity" => "planck"),
        Dict("T_of_h" => true,  "N_of_h" => true, "with_emission" => false, "albedo" => 0.0, "max_isotope_id" =>  11, "back_ground"=> 0.0, "Δλ" => 1.0e-11, "initial_intensity" => "planck"),
        Dict("T_of_h" => false, "N_of_h" => true, "with_emission" => false, "albedo" => 0.0, "max_isotope_id" =>  11, "back_ground"=> 0.0, "Δλ" => 1.0e-11, "initial_intensity" => "planck"),
    ])

function get_directory_paths()
    radtrans_root = dirname(@__DIR__)
    radtrans_root = "/home/wester/Projects/Julia/Private/SimpleRadTrans.jl"
    data_dir  = joinpath(radtrans_root, "HITRAN")
    input_dir = joinpath(radtrans_root, "radinput")
    out_root  = joinpath(radtrans_root, "radoutput")
    mkpath(input_dir)
    mkpath(joinpath(out_root, "intensity"))
    mkpath(joinpath(out_root, "spectrum"))
    data_dir, input_dir, out_root
end

function run_radition_transfer(json_file_paths)
    for json_file in json_file_paths
        spectrum = Spectrum(json_file)
        integrate(spectrum)
    end
end

function run_radition_transfer_rust(json_file_paths, rust_exe)
    for json_file in json_file_paths
        println(json_file)
        `$rust_exe $json_file`
    end
end

data_dir, input_dir, output_root = get_directory_paths()
json_file_paths = create_parameter_json(input_dir, output_root, "C", input_run_parameter)

#run_radition_transfer(json_file_paths)

rust_exe = "/home/wester/Projects/GitHub/SARM/rust/rust/target/debug/sarm"
run_radition_transfer_rust(json_file_paths, rust_exe)


