import os

import numpy as np
import os, sys
from timeit import default_timer as timer
import platform

script_root = os.path.abspath(os.path.dirname(__file__))

def interpolate(x0, y0, n):
    """Interpolate data onto a equidistant grid with n grid points

    Arguments:
        x0 {numpy float} -- [description]
        y0 {numpy float} -- [description]
        n {int} -- [description]

    Returns:
        (numpy, numpy) -- [description]
    """
    x0 = np.array(x0)
    y0 = np.array(y0)
    n0 = x0.shape[0]

    x = np.mgrid[x0[0]:x0[-1]:n*1j]
    y = np.zeros(n)

    j = 0
    for i in range(n):
        xx = x[i]
        while not (x0[j] <= xx and x[i] <= x0[j+1]):
            j += 1
            if (j > n0-2):
                break
        j = min(j, n0-2)
        v = (x[i] - x0[j]) / (x0[j+1] - x0[j])
        y[i] = (1.0 - v) * y0[j] + v * y0[j+1]

    return x, y

def make_lookup_for_Q(CO2_Q_dir, Tmax=300.0):
    """Create lokup tables for the CO2 partition function
     https://hitran.org/docs/iso-meta/
    global ID 	local ID 	Formula 	AFGL code 	Abundance 	        Molar Mass /g·mol-1 	Q(296 K) 	Q (full range) 	gi
    7 	        1 	        12C16O2 	626 	    0.984204 	        43.98983 	            286.09 	    q7.txt 	        1

    Arguments:
        CO2_Q_file {str} -- file with T, Q values (HITRAN data)
        n {int} -- number of T,Q pairs

    Returns:
        T, Q -- [description]
    """
    with open(os.path.join(CO2_Q_dir, "q7-q122-description.txt"), 'r') as ein:
        lines = ein.read().splitlines()[2:]
    paths      = []
    isotope_id = []
    isotope_c  = []
    isotope_m  = []
    gis        = []

    # mass[kg] => g/mol * mass_factor
    mass_factor = 1.0e-3/6.02214076e23

    # read the dexription file
    for i, line in enumerate(lines):
        ls   = line.split()
        global_id = int(ls[0])
        isotope_id.append(int(ls[1]))
        isotope_c.append(float(ls[4]))
        paths.append(ls[7])
        gis.append(int(ls[8]))

        mass = float(ls[5]) * mass_factor
        isotope_m.append(mass)

    T = []
    Q = []
    # read the partition function files
    for path in paths:
        with open(os.path.join(CO2_Q_dir, path), 'r') as ein:
            lines = ein.read().splitlines()
            TQ = np.array([[float(x) for x in line.split()] for line in lines])
            TT = np.array(TQ[:,0])
            dT = TT[1:]-TT[:-1]
            #if np.amax(dT) > 1.0 or np.amin(dT) < 1.0:
            #    print(TT)
            QQ = np.array(TQ[:,1])
            index = np.where(TT < Tmax, True, False)
            T.append(TT[index])
            Q.append(QQ[index])

    n = T[0].shape[0]
    m = len(T)
    TQ = np.zeros((n, m + 1), np.double)
    TQ[:,0] = T[0]
    for i in range(m):
        #if np.abs(np.sum(T[0] - T[i])) > 1.0e-3:
        #    print(np.abs(np.sum(T[0] - T[i])))
        TQ[:,i+1] = Q[i]

    return TQ, paths, isotope_id, isotope_c, isotope_m, gis


def make_T_p_over_height(data_dir):
    """Create

    Arguments:

    Returns:
        [type] -- [description]
    """

    h1 = [   0.0,   1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0,  11.0,  12.0]
    p1 = [1013.0, 902.0, 802.0, 710.0, 628.0, 554.0, 487.0, 426.0, 372.0, 324.0, 281.0, 243.0, 209.0]
    h2 = [ 13.0,  14.0,  15.0,  16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 70.0]
    p2 = [179.0, 153.0, 130.0, 111.0, 95.0, 81.2, 69.5, 59.5, 51.0, 43.4, 27.7, 13.2, 6.52, 3.33, 1.76, 0.951, 0.067]

    n    = 100
    h_p  = np.array(h1 + h2) * 1.0e3
    p_p  = np.array(p1 + p2) * 1.0e2
    h_p, p_p = interpolate(h_p, p_p, n)

    h_T = np.array([0.0, 13.0, 17.0, 25.0, 30.0, 45.0, 50.0, 70.0])*1.0e3
    T_T = np.array([288.0, 215.8, 215.7, 225.1, 233.7, 269.9, 275.7, 218.1])
    h_T, T_T = interpolate(h_T, T_T, n)

    T = np.zeros((h_T.shape[0], 2), np.double)
    p = np.zeros((h_p.shape[0], 2), np.double)

    h_T_path = os.path.join(data_dir, "h_T.npy")
    h_p_path = os.path.join(data_dir, "h_p.npy")

    T[:,0] = h_T
    T[:,1] = T_T

    p[:,0] = h_p
    p[:,1] = p_p

    np.save(h_T_path, T)
    np.save(h_p_path, p)

def make_spectrum(hitran_file, data_dir, λmin, λmax, lids, abus, masss, gis):
    h  = 6.62607004e-34
    c0 = 2.99792458e8

    with open(os.path.join(data_dir, hitran_file + ".out"), 'r') as ein:
        lines = ein.read().splitlines()[1:]
    m = np.array([[float(x) for x in line.split(',')] for line in lines])

    mid = m[:, 0]
    iid = m[:, 1]
    ν   = m[:, 2]
    ν   = ν / 1.0e-2
    λ   = 1.0 / ν

    index1 = np.where(λ >= λmin, True, False)
    index2 = np.where(λ <= λmax, True, False)
    index = index1*index2

    λ   = λ[index]
    iid = iid[index]

    S   = m[index,  3]
    A   = m[index,  4]
    γ_a = m[index,  5]
    γ_s = m[index,  6]
    ν_l = m[index,  7]
    n_a = m[index,  8]
    δ_a = m[index,  9]
    g_u = m[index, 10]
    g_l = m[index, 11]

    ν_l    = ν_l / 1.0e-2           # cm => m, bar => pascal
    γ_a    = γ_a / 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    γ_s    = γ_s / 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    δ_a    = δ_a / 1.0e-2 * 1.0e-5  # cm => m, bar => pascal
    ΔE_ul  = h * c0 / λ
    E_l    = h * c0 * ν_l
    E_u    = E_l + ΔE_ul

    n = A.shape[0]
    c = np.zeros((n, 14))

    c[:,  0] = λ
    c[:,  1] = E_l
    c[:,  2] = E_u
    c[:,  3] = S
    c[:,  4] = A
    c[:,  5] = γ_a
    c[:,  6] = γ_s
    c[:,  7] = n_a
    c[:,  8] = δ_a
    c[:,  9] = g_u
    c[:, 10] = g_l

    i = np.argmax(c[:,3])
    #print(i, c[i,:])
    # 0 1 2 3 4 5 6 7 8 9 10 11
    # 1 2 3 4 5 6 7 8 9 0 11 12

    itoj = [9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 10, 11]

    for i in range(n):
        ii = int(iid[i])
        j = itoj[ii]
        c[i, 11] = j
        c[i, 12] = masss[j]
        c[i, 13] = abus[j]

    np.save(os.path.join(data_dir, hitran_file + ".npy"), c)

    for j in range(30,50):
        s = []
        for i in range(14):
            s.append("%12.6e" % c[j, i])
        #print("   ".join(s))

def create_npy_data_files(data_dir):
    make_T_p_over_height(data_dir)

    CO2_Q_dir = os.path.join(data_dir, "CO2_Q")
    TQ, paths, isotope_id, isotope_c, isotope_m, gis = make_lookup_for_Q(CO2_Q_dir, Tmax = 300.0)
    TQ_path = os.path.join(data_dir, "T_Q.npy")
    np.save(TQ_path, TQ)

    λmin, λmax = 1.19e-5, 1.81e-5
    hitran_file = os.path.join(data_dir, "CO2_rwfmt_ISO-0-12_wl-12-18-mum")
    make_spectrum(hitran_file, data_dir, λmin, λmax, isotope_id, isotope_c, isotope_m, gis)