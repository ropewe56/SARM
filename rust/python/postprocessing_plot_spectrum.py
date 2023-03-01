import h5py as h5
import numpy as np
import pylab as plt
import sys, os
import shutil

script_root = os.path.abspath(os.path.dirname(__file__))

# https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy
# line strength
# A_ki = \dfrac{2 π e^2}{m_e c \espilon_0 \lambda^2} \dfrac{g_i}{g_k} f_{ik}
#      = \dfrac{16  \pi^3}{3 h \epsilon_0 \lambda^3 g_k} S
# \epsilon_{line} =  \dfrac{1}{4 \pi} h ν A_{ki} N_k = J/s/m³
#http://climatemodels.uchicago.edu/modtran/


def load_hdf5(input_file, group, dataset):
    """Load a hdf5 file

    Arguments:
        input_file {[type]} -- [description]
        group {[type]} -- [description]
        dataset {[type]} -- [description]

    Returns:
        [type] -- [description]
    """
    f = h5.File(input_file, "r")
    gr = f[group]
    data = gr[dataset]
    d = data[...]
    f.close()
    return d

def load(input_file, group, dataset):
    """Load hdf5 or npy file

    Arguments:
        input_file {[type]} -- [description]
        group {[type]} -- [description]
        dataset {[type]} -- [description]

    Returns:
        [type] -- [description]
    """
    if os.path.splitext(input_file)[1] == ".hdf5":
        return load_hdf5(input_file, group, dataset)
    elif os.path.splitext(input_file)[1] == ".npy":
        return np.load(input_file)

    return None

def plot_lines():
    input_file = os.path.join(script_root, "CO2.hdf5")
    group_name, dataset_name = "CO2", "data"
    data = load(input_file, group_name, dataset_name)
    # 0 wlnm
    # 1 dEnm
    # 2 Em
    # 3 Anm
    # 4 gn
    # 5 gm
    print(data[0,:])
    wl = data[:,0]
    Em = data[:,1]
    En = Em + data[:,2]
    A  = data[:,3]
    gn = data[:,4]
    L  = data[:,11]

    h = 6.62607004e-34
    c = 2.99792458e8
    kB = 1.38064852e-23
    cc = h * c * 1.0e2
    T  = 288.0
    c1 = 1.0/(4.0*np.pi) * h * c

    plt.plot(wl*1.0e6, np.exp(-En / (kB * T)) * A * gn * c1 / 2.900e2 /wl * 2.5e23 * 400.0e-6 * 1.0e3)
    plt.figure()

    input_file = os.path.join(script_root, "ea.hdf5")
    group = "wea"
    dataset = "wea"
    wea = load(input_file, group, dataset)

    plt.plot(wea[:,0]*1.0e6,wea[:,1])
    plt.figure()
    plt.plot(wea[:,0]*1.0e6,wea[:,2])

def plot_intensity_over_wavelength(ifig, path0, path, i):
    group_name, dataset_name = "spectrum", "lI"

    lI0 = load(os.path.join(path0[0], path0[1]+path0[2]), group_name, dataset_name)
    lI  = load(os.path.join(path[0],  path[1]+path[2]),   group_name, dataset_name)

    plt.figure(ifig)
    plt.plot(lI0[:,0]*1.0e6, lI0[:,1], 'r')
    plt.plot(lI[:,0]*1.0e6, lI[:,1], 'b')
    plt.xlabel("λ [μm]")
    plt.ylabel("I_λ [W/m^2/sr/m]")
    plt.savefig(os.path.join(path[0], "%d_%s_I_λ.png" % (i, path[1])))
    plt.close(ifig)

def plot_spectrum_over_wavelength(ifig1, ifig2, path, i):
    group_name, dataset_name = "conv", "conv"
    wl = load(os.path.join(path[0], path[1]+path[2]), group_name, dataset_name)

    plt.figure(ifig1)
    plt.plot(wl[:,0]*1.0e6, wl[:,1], 'b')
    plt.xlabel("λ [μm]")
    plt.ylabel("κ [a.u.]")
    plt.savefig(os.path.join(path[0], "%d_%s_conv_kap.png" % (i, path[1])))
    plt.close(ifig1)

    plt.figure(ifig2)
    plt.plot(wl[:,0]*1.0e6, wl[:,2]*1.0e-7, 'b')
    plt.xlabel("λ [μm]")
    plt.ylabel("ϵ [a.u.]")
    plt.savefig(os.path.join(path[0], "%d_%s_conv_eps.png" % (i, path[1])))
    plt.close(ifig2)


def renames(root, ext):
    spectrum, intensity, moving_average_k, planck_intensity = [], [], None, None
    dirs = set()
    for rr, dd, ff in os.walk(root):
        for f in ff:
            fs = os.path.splitext(f)
            if fs[1] == ext:
                if "ma_k" in f:
                    a = "ma_k"
                    b = "moving_average_k"
                    shutil.move(os.path.join(rr, fs[0]+fs[1]), os.path.join(rr, fs[0].replace(a, b)+fs[1]))
                elif "lI_planck" in f:
                    a = "lI_planck"
                    b = "planck_intensity"
                    shutil.move(os.path.join(rr, fs[0]+fs[1]), os.path.join(rr, fs[0].replace(a, b)+fs[1]))
                elif "conv" in f:
                    a = "conv"
                    b = "spectrum"
                    shutil.move(os.path.join(rr, fs[0]+fs[1]), os.path.join(rr, fs[0].replace(a, b)+fs[1]))
                elif "lI" in f:
                    a = "lI"
                    b = "intensity"
                    shutil.move(os.path.join(rr, fs[0]+fs[1]), os.path.join(rr, fs[0].replace(a, b)+fs[1]))


def get_all_files_by_extension(root, ext):
    spectrum, intensity, moving_average_k, planck_intensity = [], [], None, None
    dirs = set()
    for rr, dd, ff in os.walk(root):
        for f in ff:
            fs = os.path.splitext(f)
            if fs[1] == ext:
                if "moving_average_k" in f:
                    moving_average_k = (rr, fs[0], fs[1])
                elif "planck_intensity" in f:
                    planck_intensity = (rr, fs[0], fs[1])
                elif "spectrum" in f:
                    spectrum.append((rr, fs[0], fs[1]))
                    dirs.add(rr)
                elif "intensity" in f:
                    intensity.append((rr, fs[0], fs[1]))

    spectrum.sort()
    intensity.sort()
    return spectrum, intensity, moving_average_k, planck_intensity

def plot_moving_average(path):
    ma = np.load(os.path.join(path[0], path[1]+path[2]))
    print(ma.shape)
    plt.plot(ma[:,0]*1.0e6, ma[:,1])
    plt.figure()

import parameter_input

_, out_root, _, _ = parameter_input.get_root_dirs()

out_root = "/home/rolf/Projects/climate/sarm/radoutput"

subdir = 'F'
subsubdirs = os.listdir(os.path.join(out_root, subdir))
subsubdirs.sort()

#for d in subsubdirs:
#    renames(os.path.join(out_root, subdir, d), ".npy")

for d in subsubdirs:
    spectrum, intensity, moving_average_k, planck_intensity = get_all_files_by_extension(os.path.join(out_root, subdir, d), ".npy")

    if intensity != None and planck_intensity != None:
        n = len(intensity)
        for i in range(n):
            plot_spectrum_over_wavelength(1, 2, spectrum[i], i)
            plot_intensity_over_wavelength(3, planck_intensity, intensity[i], i)

