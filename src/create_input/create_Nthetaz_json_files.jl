using Common
using JSON

using PyPlot
pygui(true)
pygui(:qt5)

export get_parameter_and_create_json

"""
get_paths_dict(input_dir, out_dir)

    data_dir  -- directory where the HITRAN data are
    out_dir   -- output directory
"""
function get_paths_dict(input_dir, out_dir)
    Dict(   "out_dir"      => out_dir,
            "logfile"      => "log.log",
            "infofile"     => "info.log",
            "NCO2_path"    => joinpath(input_dir, "NCO2.npy" ),
            "theta_path"   => joinpath(input_dir, "theta_deg.npy"),
            "specdat_path" => joinpath(data_dir, "CO2_rwfmt_ISO-0-12_wl-12-18-mum.npy"),
            "T_Q_path"     => joinpath(input_dir, "T_Q.npy"  ),
            "h_p_path"     => joinpath(input_dir, "h_p.npy"  ),
            "h_T_path"     => joinpath(input_dir, "h_T.npy"  ),
            "z_path"       => joinpath(input_dir, "z.npy"    ),
            "z_iout"       => joinpath(input_dir, "z_iout.npy"),
        )
end

"""
    Get the directories where HITRAN data are,
    where to put output
    and where the rust executable is located (rust is deprecated).
"""
function get_root_dirs()
    radtrans_root = dirname(@__DIR__)
    out_root = joinpath(radtrans_root, "radoutput")
    data_dir = joinpath(radtrans_root, "HITRAN")
    exe = joinpath(radtrans_root, "radtrans_rs", "target", "release", "radtrans")
    return radtrans_root, out_root, data_dir, exe
end

"""
    Create a Vector{Float64} with z-values for the integration

    not used anymore, instead see make_z_log10
    zmin
    zmax
    dzmin
    dzmax
    n
    e

"""
function make_z_e(zmin, zmax, dzmin, dzmax, n, e)
    dz = collect(LinRange(dzmin, dzmax, n))
    dz = dz.^e
    z  = cumsum(dz)
    z  = z * zmax / maximum(z)

    z_iout = zeros(Int64, size(z,1))
    zout = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 10000.0, 70000.0]
    for zo in zout
        zz = z .- zo
        i = argmin(zz.^2)
        z[i] = zo
        z_iout[i] = 1
    end

    h = collect(LinRange(zmin, zmax, n))
    zh = exp.(h./maximum(h)*e).-1.0
    zh = zh/maximum(zh) * zmax

    z, z_iout, zh
end

"""
    make_z_log10(zmin, zmax, n)

    make a Vector of z-values: radiation transfer from z[i] to z[i+1]
    z-values are created equadistantly on a log10 scale, than z = 10^log10_z

    create a Vector{Int64} z_iout of length n, plot inetensity and spectrum where z_iou == 1

    zmin : starting z value (e.g. 0.1 m)
    zmax : end z value (70 km TAO)
    n : number of z-values

"""
function make_z_log10(zmin, zmax, n)
    log10_z = collect(LinRange(log10(max(1.0e-1, zmin)), log10(zmax), n))
    z = @. 10^log10_z

    z_iout = zeros(Int64, size(z,1))
    zout = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 7000.0, 10.0e3, 20.0e3, 30.0e3, 40.0e3, 60.0e3, 70.0e3]
    for zo in zout
        zz = z .- zo
        i = argmin(zz.^2)
        z[i] = zo
        z_iout[i] = 1
    end

    z, z_iout, z
end

"""
    Generate all parameter for the run

    input_dir : directory ehere to write input date files
    out_root : output root diectory
    subdir : the output subdirectory made from input run parameters
    input_run_parameter
"""
function wite_data_input_files(input_dir, input_run_parameter)
    par = input_run_parameter["npy_parameter"]

    # CO2, theta, npy
    write_npy(joinpath(input_dir, "NCO2.npy"),       par["NCO2"])
    write_npy(joinpath(input_dir, "theta_deg.npy"),  par["theta_deg"])

    # CO2, theta, hdf5
    save_array_as_hdf5(joinpath(input_dir, "NCO2.hdf5"),  par["NCO2"])
    save_array_as_hdf5(joinpath(input_dir, "theta_deg.hdf5"),  par["theta_deg"])

    # z, npy, hdf5
    zp = par["z"]
    #z, iout, zh = make_z_e(zp["zmin"], zp["zmax"], zp["dzmin"], zp["dzmax"], zp["n"],  zp["exponent"])
    z, iout, zh = make_z_log10(zp["zmin"], zp["zmax"], zp["n"])
    write_npy(joinpath(input_dir, "z.npy"), z)
    write_npy(joinpath(input_dir, "z_iout.npy"), iout)
    save_array_as_hdf5(joinpath(input_dir, "z.hdf5"), z)
    save_array_as_hdf5(joinpath(input_dir, "z_iout.hdf5"), iout)
end

function vector_of_input_run_parameter(input_run_parameter)
    # input_run_parameter["loop_parameter"]
    loop_parameter_vec_of_dict = []
    for (i, lp) in enumerate(input_run_parameter["loop_parameter"])
        push!(loop_parameter_vec_of_dict, lp)
    end
    return loop_parameter_vec_of_dict
end

function subsubdir_names(output_root, subdir, input_run_parameter)
    subdir_path = joinpath(output_root, subdir)
    mkpath(subdir_path)

    # input_run_parameter["loop_parameter"]
    tf = Dict(false => 0, true => 1)
    subsubdirs = []
    for (i, lp) in enumerate(input_run_parameter["loop_parameter"])
        # make subdir name
        dd = @sprintf("%02d_%d_%d_%3.1f_%d_%4.2f_%5.0e_%s", i,
                        tf[lp["T_of_h"]],
                        tf[lp["with_emission"]],
                        lp["albedo"],
                        lp["max_isotope_id"],
                        lp["back_ground"],
                        lp["Δλ"],
                        lp["initial_intensity"])
        push!(subsubdirs, (output_root, subdir, dd))
    end
    return subsubdirs
end

function create_parameter_json(input_dir, output_root, subdir, input_run_parameter)
    paths = get_paths_dict(input_dir, output_root)
    wite_data_input_files(input_dir, input_run_parameter)

    subsubdirs_vec = subsubdir_names(output_root, subdir, input_run_parameter)
    loop_parameter_vec_of_dict = vector_of_input_run_parameter(input_run_parameter)

    json_file_paths = []
    for (i, run_parameter_dict) in enumerate(loop_parameter_vec_of_dict)
        parameter_dict = Dict("run_parameter" => merge(run_parameter_dict, input_run_parameter["fixed_parameter"]))

        out_dir = joinpath(subdir, subsubdirs_vec[i][1], subsubdirs_vec[i][2], subsubdirs_vec[i][3])
        if !isdir(out_dir)
            mkpath(out_dir)
        end
        paths["out_dir"] = out_dir

        parameter_dict["paths"] = paths

        json_file_path = joinpath(out_dir, "input.json")
        open(json_file_path, "w") do out
            JSON.print(out, parameter_dict, 4)
        end
        push!(json_file_paths, json_file_path)
    end
    json_file_paths
end
