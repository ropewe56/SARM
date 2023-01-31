import numpy as np
import os, sys
from timeit import default_timer as timer
import platform

script_root = os.path.abspath(os.path.dirname(__file__))

"""
Input parameter:
    CO2 concentration values
    Inclination angle values
    Input parameter for determining height values
        Tconst = False: temperature and density values vary according to height in the atmosphere
               = True : temperature is constant, density value varies according to height in the atmosphere
        with_emission = True: include emission in the radiation transfer equation
                        False: neglect emission in the radiation transfer equation
        albedo:
        nisos: number of CO2 isotopes included
        bg: background added to the absorption and emission coefficients because of the cut off of Lorentz shape contributions
        Δλ: wavelength resolution
"""
def get_input_run_parameter():
    return {
    "NCO2"      : [400.0*1.0e-6, 800.0*1.0e-6],
    "theta_deg" : [0.0, 40.0, 80.0],
    "z"         : {"zmin": 0.0, "zmax": 7.0e4, "dzmin": 1.0, "dzmax": 1.0e3, "n": 700,  "exponent": 1},
    "loop_parameter" : [
        {"Tconst" : False, "with_emission" : True, "albedo" : 0.0, "niso" : 12, "bg": 0.03, "Δλ" : 1.0e-11},
        {"Tconst" : False, "with_emission" : True, "albedo" : 0.0, "niso" :  0, "bg": 0.03, "Δλ" : 1.0e-11},
        {"Tconst" : False, "with_emission" : False,"albedo" : 0.0, "niso" :  0, "bg": 0.03, "Δλ" : 1.0e-11},
        {"Tconst" : False, "with_emission" : True, "albedo" : 0.3, "niso" :  0, "bg": 0.03, "Δλ" : 1.0e-11},
        {"Tconst" : True,  "with_emission" : True, "albedo" : 0.0, "niso" :  0, "bg": 0.03, "Δλ" : 1.0e-11},
    ]
    }

def get_parameter_dict(data_dir):
    """Returns the dictionary holding all parameter

    Arguments:
        data_dir {str} -- HITRAN data directory

    Returns:
        dict -- parameter
    """
    h  = 6.62607004e-34
    c0 = 2.99792458e8
    kB = 1.38064852e-23
    sigma = 2.0 * np.pi**5 * kB**4  / (15.0 * h**3 * c0**2)

    parameter_dict = {
        "natconst" : {
            "h"    : h,
            "c0"   : c0,
            "kB"   : kB,
            "sigma": sigma,
            "pi"   : np.pi,
            "mCO2" : (12.011 + 2.0 * 15.999) * 1.660539040e-27
        },

        "run_parameter" :{
            "T_surface": 288.0,
            "temperature_mode"   : 1,
            "λmin"               : 1.2e-5,
            "λmax"               : 1.8e-5,
            "Δλ"                 : 1.0e-11,
            "Δλ_factor"          : 20.0,        # relative width of the line shapes
            "p_ref"              : 1013.25e2,  # N / m^2
            "T_ref"              : 296.0,      # K
            "albedo"             : 0.0,
            "with_emission"      : 1,
            "background"         : 0.0,
            "max_isotope_id"     : 10,
            "compute_F"          : True,
            "integrate"          : True,
        },
        "paths" :{
            "out_dir"      : "out_dir",
            "logfile"      : "log.log",
            "logfile_F"    : "log_F.log",
            "infofile"     : "info.log",
            "NCO2_path"    : os.path.join(data_dir, "NCO2.npy" ),
            "theta_path"   : os.path.join(data_dir, "theta_deg.npy"),
            "specdat_path"   : os.path.join(data_dir, "CO2_rwfmt_ISO-0-12_wl-12-18-mum.npy"),
            "T_Q_path"     : os.path.join(data_dir, "T_Q.npy"  ),
            "h_p_path"     : os.path.join(data_dir, "h_p.npy"  ),
            "h_T_path"     : os.path.join(data_dir, "h_T.npy"  ),
            "z_path"       : os.path.join(data_dir, "z.npy"    ),
            "z_iout"       : os.path.join(data_dir, "z_iout.npy"    ),
        }
    }
    return parameter_dict

def get_root_dirs():
    """Get the directories where HITRAN data are, where to put output and where the rust executable is located.

    Returns:
        tuple -- [description]
    """
    sarm_root = os.path.abspath(os.path.join(script_root, ".."))
    exe = os.path.join(sarm_root, "sarm_rs", "target", "release", "sarm")
    return sarm_root, exe

def make_z(zmin, zmax, dzmin, dzmax, n, e):
    """Create a numpy array with z-values for the integration

    Arguments:
        zmin {float} -- [description]
        zmax {float} -- [description]
        dzmin {float} -- [description]
        dzmax {[type]} -- [description]
        n {float} -- [description]
        e {float} -- [description]

    Returns:
        [type] -- [description]
    """
    dz = np.mgrid[dzmin:dzmax:1j*n]
    dz = dz**e
    z  = np.cumsum(dz)
    z  = z * zmax / np.amax(z)

    z_iout = np.zeros(z.shape[0], np.int32)
    zout = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 10000.0, 70000.0]
    for zo in zout:
        zz = z - zo
        i = np.argmin(zz**2)
        z[i] = zo
        z_iout[i] = 1

    return z, z_iout

def get_all_parameter(data_dir, output_root, subdir):
    """Generate all parameter for the run

    Arguments:
        subdir {[type]} -- [description]

    Returns:
        [type] -- [description]
    """


    # paths to directories and the rust_exe
    sarm_root, rust_exe = get_root_dirs()

    parameter_dict = get_parameter_dict(data_dir)

    subdir_path = os.path.join(output_root, subdir)
    if not os.path.isdir(subdir_path):
        os.makedirs(subdir_path)

    input_run_parameter = get_input_run_parameter()

    np.save(os.path.join(data_dir, "NCO2.npy"),      np.array(input_run_parameter["NCO2"]))
    np.save(os.path.join(data_dir, "theta_deg.npy"), np.array(input_run_parameter["theta_deg"]))

    zp = input_run_parameter["z"]
    z, iout = make_z(zp["zmin"], zp["zmax"], zp["dzmin"], zp["dzmax"], zp["n"],  zp["exponent"])
    np.save(os.path.join(data_dir, "z.npy"), z)
    np.save(os.path.join(data_dir, "z_iout.npy"), iout)

    run_parameter = []
    dirs = []
    for i, lp in enumerate(input_run_parameter["loop_parameter"]):

        #{"Tconst" : False, "Nconst": False, "with_emission" : True, "intensity" : 1.0, "niso" : 12, "bg": 0.03, "Δλ" : 1.0e-11},

        dd = "%02d_%d_%d_%3.1f_%d_%4.2f_%5.0e" % (i, int(lp["Tconst"]), int(lp["with_emission"]),
                                                     lp["albedo"], lp["niso"], lp["bg"],
                                                     lp["Δλ"])

        run_parameter.append([int(lp["Tconst"]), int(lp["with_emission"]),
                                  lp["albedo"], lp["niso"], lp["bg"], lp["Δλ"]])

        dirs.append([output_root, subdir, dd])

    return parameter_dict, run_parameter, dirs, rust_exe, data_dir

