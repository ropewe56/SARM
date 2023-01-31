import numpy as np
import os, sys
import pylab as plt
import json
import subprocess
import platform
import logging
import threading
from timeit import default_timer as timer

import parameter_input
import create_npy_files

format = "%(asctime)s: %(message)s"
logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")
script_root = os.path.abspath(os.path.dirname(__file__))

ext_root = os.path.abspath(os.path.join(script_root, "..", "rust_py", "target", "release"))
sys.path.insert(0, ext_root)
if not os.path.exists(os.path.join(ext_root, "pysarm.so")) and \
       os.path.exists(os.path.join(ext_root, "libpysarm.so")):
    import shutil
    shutil.move(os.path.exists(os.path.join(ext_root, "libpysarm.so")),
                os.path.exists(os.path.join(ext_root, "pysarm.so")))


def create_json(parameter, par, dirs):
    """Create json from parameter dict

    Arguments:
        parameter {dict} -- parameter dictionary
        par {list} -- list of run parameter
        dirs {list of str} -- [description]
        exe {str} -- [description]
    """
    out_dir = os.path.join(dirs[0], dirs[1], dirs[2])
    parameter["paths"]["out_dir"]   = out_dir

    parameter["run_parameter"]["temperature_mode"] = par[0]
    parameter["run_parameter"]["with_emission"]    = par[1]
    parameter["run_parameter"]["albedo"]           = par[2]
    parameter["run_parameter"]["max_isotope_id"]   = par[3]
    parameter["run_parameter"]["background"]       = par[4]
    parameter["run_parameter"]["Δλ"]               = par[5]

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    jsonfile = os.path.join(out_dir, 'input.json')
    with open(jsonfile, 'w') as out:
        out.write(json.dumps(parameter, sort_keys=True, indent=4))
    #print("jsonfile = ", jsonfile)
    return jsonfile

def run_rust_exe(jsonfile, exe):
    arg = "--json_file=%s" % jsonfile
    print("exe = %s\narg = %s" % (exe, arg))
    subprocess.run([exe, arg])


def run():
    data_dir = "/home/wester/Projects/Julia/Private/SimpleRadTrans.jl/radinput"
    parameter_dict, run_parameter, dirs, rust_exe, data_dir = parameter_input.get_all_parameter(data_dir, "E")

    create_npy_files.create_npy_data_files(data_dir)

    json_files = []
    for i,par in enumerate(run_parameter):
        #print(par, dirs[i])
        json_file = create_json(parameter_dict, par, dirs[i])
        json_files.append(json_file)

    import pysarm
    for json_file in json_files:
        pysarm.call_spectrum(json_file)
        #run_rust_exe(json_file, rust_exe)
        print("======================================================================")

start = timer()
run()
end = timer()
print("time = %5.2f min" % ((end - start)/60.0))
