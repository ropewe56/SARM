use ndarray::prelude::*;
use crate::spectrum::ioutil::{save_npy_2, save_npy_1, load_npy_2, load_npy_1, open_file};

pub struct Line {
    pub λ_ul0      : f64,
    pub E_l        : f64,
    pub E_u        : f64,
    pub S          : f64,
    pub A_ul       : f64,
    pub γ_a        : f64,
    pub γ_s        : f64,
    pub n_a        : f64,
    pub δ_a        : f64,
    pub g_u        : f64,
    pub g_l        : f64,

    pub λ_ul       : f64,
    pub ΔλL        : f64,
    pub ΔλL0       : f64,
    pub ΔλG        : f64,
    pub B_ul       : f64,
    pub B_lu       : f64,
    pub f_u        : f64,
    pub κ          : f64,
    pub ϵ          : f64,
    pub N_l        : f64,
    pub N_u        : f64,

    pub κ_sum      : f64,

    pub isotope_id : usize,
    pub isotope_c  : f64,
    pub isotope_m  : f64,
}

pub fn save_lines_data(lines: &Vec<Line>, fname: &str)  {
    let nb_lines = lines.len();
    let mut a : Array2<f64> = Array2::<f64>::zeros((nb_lines, 17));
    for iline in 0..nb_lines {
        a[[iline, 0]] = lines[iline].λ_ul0   ;
        a[[iline, 1]] = lines[iline].E_l     ;
        a[[iline, 2]] = lines[iline].E_u     ;
        a[[iline, 3]] = lines[iline].S       ;
        a[[iline, 4]] = lines[iline].A_ul    ;
        a[[iline, 5]] = lines[iline].γ_a     ;
        a[[iline, 6]] = lines[iline].γ_s     ;
        a[[iline, 7]] = lines[iline].n_a     ;
        a[[iline, 8]] = lines[iline].δ_a     ;
        a[[iline, 9]] = lines[iline].g_u     ;
        a[[iline,10]] = lines[iline].g_l     ;

        a[[iline,11]] = lines[iline].λ_ul    ;
        a[[iline,12]] = lines[iline].ΔλL     ;
        a[[iline,13]] = lines[iline].ΔλG     ;
        a[[iline,14]] = lines[iline].B_ul    ;
        a[[iline,15]] = lines[iline].B_lu    ;
        a[[iline,16]] = lines[iline].κ       ;
        a[[iline,17]] = lines[iline].ϵ       ;
        a[[iline,18]] = lines[iline].N_l     ;
        a[[iline,19]] = lines[iline].N_u     ;
    }

    save_npy_2(fname, &a);
    //crate::util::hdf5_util::save_hdf5_2(fname, &a, "lines", "lines");
}
