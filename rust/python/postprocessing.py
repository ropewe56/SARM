import os
import numpy as np
import pylab as plt

import parameter_input

def compute_intensity(z, I, κ, ϵ):
    """[summary]
    
    Arguments:
        z {array} -- [description]
        I {array} -- [description]
        κ {array} -- [description]
        ϵ {array} -- [description]
    Returns:
        {array} -- intensity
    """
    dz = 100.0
    exp_κ = np.exp(-κ * dz)
    a     = ϵ / κ * (1.0 - exp_κ)
    b     = ϵ * dz
    I_    = np.zeros(z.shape[0])
    I_[0] = I[0]
    limit = 0.01
    ee   = 0.0
    for i in range(z.shape[0]-1):
        ee = ee + ϵ[i] * dz
        if κ[i] * dz < limit:
            I_[i+1] = I_[i] * exp_κ[i] + ϵ[i] * dz
        else:
            I_[i+1] = I_[i] * exp_κ[i] + a[i]
    #print(ee)
    return I_


class DataSets:
    """Read all result data from the log.log files
    
    Returns:
        [type] -- [description]
    """
    def __init__(self, out_root, subdir, subsubdirs):
        self.dirs = []
        for ssd in subsubdirs:
            print(ssd)
            self.dirs.append(os.path.join(out_root, subdir, ssd))

        self.dirs.sort()

        self.datasets = {}

        for id, dd in enumerate(self.dirs):
            self.read_data(id, dd)

    def get_size(self):
        return len(self.datasets )

    def read_data(self, id, result_dir):
        """[summary]
        
        Arguments:
            id {[type]} -- [description]
            result_dir {[type]} -- [description]
        """
        with open(os.path.join(result_dir, "log.log")) as ein:
            lines = ein.read().splitlines()
        ls = lines[0].strip('#').split()
        print(len(lines), result_dir)

        nN   = int(ls[0])     # nb NCo2
        nθ   = int(ls[1])     # nb theta
        nz   = int(ls[2])+1   # nb z, there are two iz = 0 rows in log.log file
        ncol = int(ls[3])     # nb columns
        #        0     1  2  3  4  5  6  7  8    9 
        #    iz, z, NCO2, θ, I, T, N, ϵ, κ, ΔλL, ΔλD

        # 3  2  701 10
        M = np.zeros((nN, nθ, nz, ncol))
        for i,line in enumerate(lines[1:]):
            if line[0] == '#':
                if i > 0:
                    M[iN, iθ, :, :] = np.array(m)
                c = line[1:].split()
                iN = int(c[0])
                iθ = int(c[1])
                m = []
            else:
                ls = line.split(',')
                a = []
                for entry in ls[1:]:
                    e = entry.split("=")
                    a.append(float(e[1]))
                m.append(a)
        #print(len(m), M.shape)
        M[iN, iθ, :, :] = np.array(m)
        self.datasets[id] = M

    def get_last(self):
        return self.datasets[-1]

    def get(self, id, iN, iθ, ie):
        """Get data 
        
        Arguments:
            id {int} -- subsubdir
            iN {int} -- NCO2
            iθ {int} -- theta
            ie {int} -- column
                0     1  2  3  4  5  6  7  8    9 
            iz, z, NCO2, θ, I, T, N, ϵ, κ, ΔλL, ΔλD
        
        Returns:
            array -- a column
        """
        ds = self.datasets[id]
        #       NCO2, theta, z, icol
        return ds[iN, iθ, :, ie]

# indices
iz = 0
iC = 1
iθ = 2
iI = 3
iT = 4
iN = 5
iϵ = 6
iκ = 7

def integrate_intensity(data_sets, id, nθ, iN, NCO2, color1, color2):
    """Integrate intensity ove angle theta
    
    Arguments:
        data_sets {[type]} -- [description]
        id {[type]} -- [description]
        nθ {[type]} -- [description]
        iN {[type]} -- [description]
        NCO2 {[type]} -- [description]
        color1 {[type]} -- [description]
        color2 {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    import scipy.integrate
    θ_0 = np.deg2rad(data_sets.get(id, iN, 0, iθ)[0]) # theta
    θ_1 = np.deg2rad(data_sets.get(id, iN, 1, iθ)[0]) # theta
    θ_2 = np.deg2rad(data_sets.get(id, iN, 2, iθ)[0]) # theta
    I_0 =            data_sets.get(id, iN, 0, iI)[-1] # intensity at TOA
    I_1 =            data_sets.get(id, iN, 1, iI)[-1] # intensity at TOA
    I_2 =            data_sets.get(id, iN, 2, iI)[-1] # intensity at TOA

    # qubic approximation of I(θ)
    R1 = I_1 - I_0
    R2 = I_2 - I_0

    a0 = I_0
    det = θ_1**2 * θ_2**3 - θ_1**3 * θ_2**2
    a2 = (R1 * θ_2**3 - R2 * θ_1**3) / det
    a3 = (R2 * θ_1**2 - R1 * θ_2**2) / det

    c1 = scipy.integrate.quad(lambda x: np.cos(x)*np.sin(x),      0.0, np.pi*0.5)# θ_2)
    c2 = scipy.integrate.quad(lambda x: np.cos(x)*np.sin(x)*x**2, 0.0, np.pi*0.5)# θ_2)
    c3 = scipy.integrate.quad(lambda x: np.cos(x)*np.sin(x)*x**3, 0.0, np.pi*0.5)# θ_2)

    θ = np.mgrid[0.0:np.pi*0.5:100j]
    I = a0 + a2*θ**2 + a3*θ**3

    # plot
    plt.plot(θ, I, color1, label='%d ppm, cubic approximation' % NCO2)
    plt.plot([θ_0, θ_1, θ_2], [I_0, I_1, I_2], color2+'o', label="%d ppm, computed" % NCO2)
    plt.xlabel("angle θ [rad]")
    plt.ylabel("TOA flux I(θ) [W/m²]")
    plt.legend(loc='best')

    # integrated intensity
    Iint = 2.0*np.pi * (a0*c1[0] + a2*c2[0] + a3*c3[0])

    return Iint

def plot_results(data_sets, id):

    plotfunc = plt.semilogx
    tickfunc = plt.xticks

    plotfunc = plt.plot
    tickfunc = lambda x: x

    root = data_sets.dirs[id]

    z = data_sets.get(id, 0, 0, iz)

    iN, nθ = 0, 2
    plt.figure(1)
    II_0 = integrate_intensity(data_sets, id, nθ, iN,   400, 'b', 'r')
    II_1 = integrate_intensity(data_sets, id, nθ, iN+1, 800, 'g', 'm')
    plt.savefig(os.path.join(root, "cubic_approximation.png"))
    plt.close(1)

    with open(os.path.join(root, "result.txt"), 'w') as out:
        out.write("I_400 - I_800 = %12.6e\n" % (II_0 - II_1))
    print("I_400 - I_800 = %12.6e, %12.6e, %12.6e" % (II_0 - II_1, II_0, II_1))

    for j, zz in enumerate(z):
        if zz > 0.2:
            break
    j = 1
    plt.figure(2)
    plotfunc(z[j:], data_sets.get(id, 0, 0, iI)[j:],  'b', label="400 ppm")
    plotfunc(z[j:], data_sets.get(id, 1, 0, iI)[j:],  'g', label="800 ppm")
    plt.xlabel("h [m]")
    plt.ylabel("I [W/m²]")
    plt.legend(loc='best')
    tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #plt.title(dd[2])
    plt.savefig(os.path.join(root, "I_vs_h.png"))
    #0 : 10.0, 1: 6.4, 2: 10.0, 3: 10, 4: 10, 5: 6.37, 6: 10, 9: 9.12
    plt.close(2)

    # ϵ Tv
    plt.figure(3)
    plotfunc(z[j:], ((data_sets.get(id, 0, 0, iϵ) - data_sets.get(id, 0, 0, iκ)))[j:]*1.0e3,  'b', label="400 ppm")
    plotfunc(z[j:], ((data_sets.get(id, 1, 0, iϵ) - data_sets.get(id, 1, 0, iκ)))[j:]*1.0e3,  'g', label="800 ppm")
    plt.xlabel("h [m]")
    plt.ylabel("(ϵ - κI) [mW/m³]")
    tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #plt.title(dd[2])
    plt.legend(loc='best')
    plt.savefig(os.path.join(root, "ϵ-κI_vs_h.png"))
    plt.close(3)

    plt.figure(4)
    plotfunc(z[j:], data_sets.get(id, 0, 0, iϵ)[j:],  'b', label="400 ppm")
    plotfunc(z[j:], data_sets.get(id, 1, 0, iϵ)[j:],  'g', label="800 ppm")
    plt.xlabel("h [m]")
    plt.ylabel("ϵ [W/m³]")
    tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #plt.title(dd[2])
    plt.legend(loc='best')
    plt.savefig(os.path.join(root, "ϵ_vs_h.png"))
    plt.close(4)

    plt.figure(5)
    plotfunc(z[j:], data_sets.get(id, 0, 0, iκ)[j:],  'b', label="400 ppm")
    plotfunc(z[j:], data_sets.get(id, 1, 0, iκ)[j:],  'g', label="800 ppm")
    plt.xlabel("h [m]")
    plt.ylabel("κI [W/m³]")
    tickfunc([1.0, 10.0, 1.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 7.0e4])
    #plt.title(dd[2])
    plt.legend(loc='best')
    plt.savefig(os.path.join(root, "κI_vs_h.png"))
    plt.close(6)

_, out_root, _, _ = parameter_input.get_root_dirs()

out_root = "/home/rolf/Projects/climate/sarm/radoutput"

subdir = "F"
subsubdirs = os.listdir(os.path.join(out_root, subdir))
subsubdirs.sort()
subsubdirs = subsubdirs[0:4]

data_sets = DataSets(out_root, subdir, subsubdirs)

for id in data_sets.datasets.keys():
    plot_results(data_sets, id)

#plt.show()

