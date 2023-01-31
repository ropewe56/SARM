import os

import numpy as np
import os, sys
from timeit import default_timer as timer
import platform

script_root = os.path.abspath(os.path.dirname(__file__))

import parameter_input


def make_z(zmin, zmax, dzmin, dzmax, n, e, c):
    """Create a numpy array with z-values for the integration

    Arguments:
        zmin {float} -- [description]
        zmax {float} -- [description]
        dzmin {float} -- [description]
        dzmax {[type]} -- [description]
        n {float} -- [description]
        e {float} -- [description]
        c {float} -- [description]

    Returns:
        [type] -- [description]
    """
    dz = np.mgrid[dzmin:dzmax:1j*n]
    dz = dz**e
    z  = np.cumsum(dz)
    z  = z * zmax / np.amax(z)
    if c != None:
        dz[0:-1] = z[1:] - z[0:-1]
        dz[-1] = dz[-2]
    return z

def get_all_parameter():
    """Generate all parameter for the run

    Arguments:
        subdir {[type]} -- [description]

    Returns:
        [type] -- [description]
    """

    input_run_parameter = parameter_input.input_run_parameter

    # paths to directories and the rust_exe
    sarm_root, out_root, data_dir, rust_exe = parameter_input.get_root_dirs()

    parameter_dict = parameter_input.get_parameter_dict(data_dir)

    subdir_path = os.path.join(out_root, input_run_parameter["subdir"])
    if not os.path.isdir(subdir_path):
        os.makedirs(subdir_path)

    np.save(os.path.join(data_dir, "NCO2.npy"),      np.array(input_run_parameter["NCO2"]))
    np.save(os.path.join(data_dir, "theta_deg.npy"), np.array(input_run_parameter["theta_deg"]))

    zp = input_run_parameter["z"]
    z = make_z(zp["zmin"], zp["zmax"], zp["dzmin"], zp["dzmin"], zp["n"],  zp["exponent"], zp["c"])
    np.save(os.path.join(data_dir, "z.npy"), z)

    run_parameter = []
    dirs = []
    for i, lp in enumerate(input_run_parameter["loop_parameter"]):

        #{"Tconst" : False, "Nconst": False, "with_emission" : True, "intensity" : 1.0, "niso" : 12, "bg": 0.03, "Δλ" : 1.0e-11},

        dd = "%02d_%d_%d_%3.1f_%d_%4.2f_%5.0e" % (i, int(lp["Tconst"]), int(lp["Nconst"]), int(lp["with_emission"]),
                                                         lp["intensity"], lp["niso"], lp["bg"], lp["Δλ"])

        run_parameter.append([int(lp["Tconst"]), int(lp["Nconst"]), int(lp["with_emission"]),
                                  lp["intensity"], lp["niso"], lp["bg"], lp["Δλ"]])

        dirs.append([out_root, input_run_parameter["subdir"], dd])

    return parameter_dict, run_parameter, dirs, rust_exe, data_dir

