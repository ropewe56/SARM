
use ndarray::prelude::*;
//use rayon::prelude::*;
//use ndarray_parallel::prelude::*;
//use ndarray::Zip;
use ndarray::{Data, RemoveAxis, Zip};
use std::f64::consts::PI;

use std::cmp;
use std::time::{Duration, Instant};
use std::fmt::format;

use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
//use std::path::Path;
//use std::path::PathBuf;

use std::collections::HashMap;

use crate::spectrum::resultdata::ResultData;
use crate::spectrum::parameter::Parameter;
use crate::spectrum::profiles::VoigtProfile;
use crate::spectrum::ioutil::{save_npy_2, save_npy_1, load_npy_2, load_npy_1, load_npy_1_i, open_file};
use crate::spectrum::rdutil::{float_mean, float_max, moving_average};
use crate::spectrum::spectralline::{Line};
use crate::spectrum::interpolator::{Interpolator};

/// The Spectrum struct
pub struct Spectrum {
    /// reference temperature
    pub T_ref    : f64,
    /// reference pressure
    pub p_ref    : f64,
    /// vector of line data
    pub lines    : Vec<Line>,
    /// the Interplator for the
    pub Q_CO2    : Interpolator,
    /// pressure over height
    pub p_vs_h   : Interpolator,
    /// temperature over height
    pub t_vs_h   : Interpolator,
    /// CO2 concentration values to iterate over
    pub NCO2     : Array1<f64>,
    /// the angles to iterate over
    pub theta_deg: Array1<f64>,
    /// height values
    pub z        : Array1<f64>,
    /// whwre to save intensites
    pub z_iout   : Array1<i64>,
    /// wavelengths
    pub λ        : Array1<f64>,
    /// absorptions over wavelengths
    pub κ_c      : Array1<f64>,
    /// emissivity  over wavelengths
    pub ϵ_c      : Array1<f64>,
    /// Parameter object
    pub par      : Parameter,

    /// logfile
    pub logfile  : File,
    /// logfile F
    pub logfile_F: File,
    pub resultdata:   ResultData,
}

/// Spectrum implementation
impl Spectrum {

    /// Spectrum constructor
    /// json_file is a json string containing all necessary information
    pub fn new(json_file: &str) -> Spectrum {
        // create a Parameter object form the jso strinh
        let par : Parameter = Parameter::new(json_file);
        par.print_parameter();

        // spectral data (HITRAN)
        let spectral_data = load_npy_2(&par.specdat_path);
        // height values
        let z             = load_npy_1(&par.z_path);
        let z_iout        = load_npy_1_i(&par.z_iout);
        // CO2 concentration values
        let NCO2          = load_npy_1(&par.NCO2_path);
        // angles
        let theta_deg     = load_npy_1(&par.theta_path);

        // CO2 partition functions
        let Q_CO2  : Interpolator = Interpolator::new(&par.T_Q_path);
        // pressure over height
        let p_vs_h : Interpolator = Interpolator::new(&par.h_p_path);
        // temperature over height
        let t_vs_h : Interpolator = Interpolator::new(&par.h_T_path);

        // spectral data
        // number of lines
        let nb_lines = spectral_data.nrows();
        // number of data per line
        let nb_par = spectral_data.ncols();

        // create Vec of nb_lines
        let mut lines : Vec<Line> = Vec::with_capacity(nb_lines);

        let mut isotope_ids = HashMap::new();

        // Create Line objects from iput line data
        for iline in 0..nb_lines {
            let λ_ul0  = spectral_data[[iline, 0]];
            let E_l    = spectral_data[[iline, 1]];
            let E_u    = spectral_data[[iline, 2]];
            let S      = spectral_data[[iline, 3]];
            let A_ul   = spectral_data[[iline, 4]];
            let γ_a    = spectral_data[[iline, 5]];
            let γ_s    = spectral_data[[iline, 6]];
            let n_a    = spectral_data[[iline, 7]];
            let δ_a    = spectral_data[[iline, 8]];
            let g_u    = spectral_data[[iline, 9]];
            let g_l    = spectral_data[[iline,10]];
            let iso_id = spectral_data[[iline,11]];
            let iso_m  = spectral_data[[iline,12]];
            let iso_c  = spectral_data[[iline,13]];

            // Einstein coefficient of induced emission
            let B_ul = A_ul * λ_ul0 * λ_ul0 * λ_ul0 / (8.0*PI * par.h);
            // Einstein coefficient of absorption
            let B_lu = g_u / g_l * B_ul;

            // Line pressure broadening coeffcient
            let γ = (γ_a * 1.0e5);
            // Lorentz line width
            let ΔλL0 = λ_ul0*λ_ul0 * γ;

            // id of
            let isotope_id = iso_id as usize;
            let h = isotope_ids.entry(isotope_id).or_insert(0);
            *h += 1;

            // Add Line only if isotope id <= max isotope id
            if isotope_id <= par.max_isotope_id {
                lines.push(Line{
                    λ_ul0      : λ_ul0 ,
                    E_l        : E_l ,
                    E_u        : E_u ,
                    S          : S   ,
                    A_ul       : A_ul,
                    γ_a        : γ_a ,
                    γ_s        : γ_s ,
                    n_a        : n_a ,
                    δ_a        : δ_a ,
                    g_u        : g_u ,
                    g_l        : g_l ,

                    λ_ul       : λ_ul0,
                    ΔλL        : ΔλL0,
                    ΔλL0       : ΔλL0,
                    ΔλG        : 0.0,

                    B_ul       : B_ul,
                    B_lu       : B_lu,
                    f_u        : ΔλL0,
                    κ_sum      : 0.0,

                    κ          : 0.0,
                    ϵ          : 0.0,

                    N_l        : 0.0,
                    N_u        : 0.0,

                    isotope_id : isotope_id,
                    isotope_c  : iso_c,
                    isotope_m  : iso_m
                } );
            }
        }

        // wavelengths array
        let nb_λ = ((par.λmax - par.λmin) / par.Δλ) as usize;
        let λ : Array1<f64> = Array::linspace(par.λmin, par.λmax, nb_λ);

        // logfile
        let mut logfile = open_file(par.out_dir.as_str(), par.logfile.as_str());
        let mut logfile_F = open_file(par.out_dir.as_str(), par.logfile_F.as_str());

        //resukts
        let mut resultdata = ResultData::new();

        // output of isotope ids
        let mut keys_vec: Vec<(&usize, &usize)> = isotope_ids.iter().collect();
        keys_vec.sort_by(|a, b| a.0.cmp(b.0));

        // infofile
        let mut infofile = open_file(par.out_dir.as_str(), par.infofile.as_str());

        let mut out = Vec::new();
        for (id, nb) in &keys_vec {
            out.push(format!("id = {:3}, nb = {:5}\n", id, nb));
        }
        out.push(format!("nb_lines used = {:3}, par.max_isotope_id = {}\n", lines.len(), par.max_isotope_id));
        out.push(format!("background    = {:5.3}\n", par.background));
        let joined = out.join("");
        match infofile.write_all(joined.as_bytes()) {
            Err(why) => { panic!("couldn't write to infofile: {}", why.to_string()) },
            Ok(_) => 0
        };
        print!("{}", joined);
        infofile.flush();

        Spectrum{T_ref    : par.T_ref,
                 p_ref    : par.p_ref,
                 lines    : lines,
                 Q_CO2    : Q_CO2,
                 p_vs_h   : p_vs_h,
                 t_vs_h   : t_vs_h,
                 NCO2     : NCO2,
                 theta_deg: theta_deg,
                 par      : par,
                 z        : z,
                 z_iout   : z_iout,
                 λ        : λ,
                 κ_c      : Array1::<f64>::zeros(nb_λ),
                 ϵ_c      : Array1::<f64>::zeros(nb_λ),
                 logfile  : logfile,
                 logfile_F: logfile_F,
                 resultdata  : resultdata,
            }
    }

    /// Save convolved absorption and emission spectra to hdf5 file
    pub fn save_convolved_to_hdf5(&self, path: &str, group_name: &str, dataset_name: &str)  {
        let nb_λ = self.λ.len();
        let mut a : Array2<f64> = Array2::<f64>::zeros((nb_λ, 3));
        for i in 0..nb_λ {
            a[[i, 0]] = self.λ[i];
            a[[i, 1]] = self.κ_c[i];
            a[[i, 2]] = self.ϵ_c[i];
        }
        let fname = [&path, ".hdf5"].join("");
        //crate::util::hdf5_util::save_hdf5_2(fname, &a, &gr_name, &ds_name);
    }

    /// Save convolved absorption and emission spectra to ny file
    pub fn save_spectrum_as_npy(&self, path: &str)  {
        let nb_λ = self.λ.len();
        let mut a : Array2<f64> = Array2::<f64>::zeros((nb_λ, 3));
        for i in 0..nb_λ {
            a[[i, 0]] = self.λ[i];
            a[[i, 1]] = self.κ_c[i];
            a[[i, 2]] = self.ϵ_c[i];
        }
        let fname = [&path, ".npy"].join("");
        save_npy_2(&fname, &a);
    }

    /// Save wavelength and intensity to npy file
    ///
    /// λ - wavelengths array
    ///
    /// I - intensity array
    ///
    /// path - filepath
    pub fn save_intensity_as_npy(&self, λ: &Array1<f64>, I: &Array1<f64>, path: &str) {
        let  n = I.len();
        let mut A : Array2<f64> = Array2::<f64>::zeros((n, 2));
        for i in 0..n {
            A[[i,0]] = λ[i];
            A[[i,1]] = I[i];
        }
        let fname = [&path, ".npy"].join("");
        save_npy_2(&fname, &A);
    }

    pub fn compute_and_save_planck(&self, λmin: f64, λmax: f64, nl: usize, T: f64, hdf5_path: &str) {
        let λ0 : Array1<f64> = Array::linspace(λmin, λmax, nl);
        let IPlanck = self.planck_λ(T, &λ0);
        self.save_intensity_as_npy(&λ0, &IPlanck, &hdf5_path);
    }

    // Planck intensity at temperature temp and wavelength λ
    pub fn planck(&self, temp: f64, λ: f64) -> f64 {
        // I = 2.0 * π * h * c**2 / (λ**5) / ( exp( (h * c / (kB * T *λ) ) - 1.0 )
        let a = 2.0 * PI * self.par.h * self.par.c0*self.par.c0;
        let b = self.par.h * self.par.c0 / temp / self.par.kB;
        let planck = a / (λ*λ*λ*λ*λ) / ((b/λ).exp() - 1.0);
        planck
    }

    // Planck intensities at temperature temp and wavelengths λ
    pub fn planck_λ(&self, temp: f64, λ: &Array1<f64>) -> Array1<f64> {
        // I = 2.0 * π * h * c**2 / (λ**5) / ( exp( (h * c / (kB * T *λ) ) - 1.0 )
        let a = 2.0 * PI * self.par.h * self.par.c0*self.par.c0;
        let b = self.par.h * self.par.c0 / temp / self.par.kB;
        let planck = λ.mapv(|x| a / (x*x*x*x*x) / ((b/x).exp() - 1.0));
        planck
    }

    /// Compute emission and absorption coefficients of the lines
    ///
    /// T - temperature
    ///
    /// N - atmosphere density
    ///
    /// p - pressure
    ///
    /// NCO2 - CO2 concentration
    pub fn compute_emission_and_absorption(&mut self,
                    T: f64, N: f64, p: f64, NCO2: f64) ->  (f64, f64) {

        let mut flag : bool = false;

        let n1 = self.lines.len();

        let dΩ = 1.0;

        let β = 1.0/(self.par.kB * T);

        let mut ΔλL_mean : f64 = 0.0;
        let mut ΔλD_mean : f64 = 0.0;

        for iline in 0..n1 {
            let λ_ul0  = self.lines[iline].λ_ul;
            let n_a    = self.lines[iline].n_a;
            let δ_a    = self.lines[iline].δ_a;
            let γ_a    = self.lines[iline].γ_a;
            let γ_s    = self.lines[iline].γ_s;
            let E_l    = self.lines[iline].E_l ;
            let E_u    = self.lines[iline].E_u ;
            let A_ul   = self.lines[iline].A_ul;
            let g_l    = self.lines[iline].g_l ;
            let g_u    = self.lines[iline].g_u ;
            let iso_id = self.lines[iline].isotope_id;
            let iso_c  = self.lines[iline].isotope_c ;
            let iso_m  = self.lines[iline].isotope_m ;

            let Q     = self.Q_CO2.yy_at(T, iso_id);
            let Niso  = N * NCO2 * iso_c;

            let λ_ul = λ_ul0 / (1.0 + λ_ul0 * δ_a * p);
            self.lines[iline].λ_ul = λ_ul;


            // γ = (Tref/T)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})
            let dT = f64::powf(self.par.T_ref/T, n_a);
            let γ = dT * (γ_a * p * (1.0 - NCO2) + γ_s * p * NCO2);

            self.lines[iline].ΔλL = λ_ul*λ_ul * γ;
            self.lines[iline].ΔλG = (2.0 * self.par.kB * T / iso_m).sqrt() / self.par.c0 * λ_ul;
            //println!("{:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e}", γ, γ_s, γ_a, dT, self.lines[iline].ΔλL, self.lines[iline].ΔλG, λ_ul - λ_ul0);

            ΔλL_mean += self.lines[iline].ΔλL;
            ΔλD_mean += self.lines[iline].ΔλG;

            //println!("{:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e}", γ, γ_a, dT, self.lines[iline].Δλ, λ_ul - λ_ul0);


            let B_ul = A_ul * λ_ul * λ_ul * λ_ul / (8.0*PI * self.par.h);
            let B_lu = g_u / g_l * B_ul;

            let N_l  = g_l * (- E_l * β).exp() / Q * Niso;
            let N_u  = g_u * (- E_u * β).exp() / Q * Niso;
            self.lines[iline].N_l = N_l;
            self.lines[iline].N_u = N_u;

            // ϵ_λ = h * c / λ_0 / (4 * π) * N_u * A_ul * f_λ
            self.lines[iline].ϵ = self.par.h * self.par.c0 / λ_ul *  N_u * A_ul * dΩ / (4.0 * PI);

            // κ_λ = h / λ_0 * N_l * B_lu * (1 - N_u/N_l * g_l/g_u) * λ_0**2 / c * f_λ
            self.lines[iline].κ = self.par.h * λ_ul / self.par.c0 * B_lu * N_l - self.par.h * λ_ul / self.par.c0 * B_ul * N_u;
        }
        // return average line widths
        (ΔλL_mean / (n1 as f64), ΔλD_mean / (n1 as f64))
    }

    /// Convolve the all lines using their line shape
    ///
    /// T - temperature
    ///
    /// N - density
    pub fn convolution(&mut self, T: f64, N: f64)  {
        let nb_lines = self.lines.len();
        let nb_λ    = self.λ.len();

        let λ_1 : f64 = self.λ[0];
        let λ_2 : f64 = self.λ[nb_λ - 1];
        let Δλ  : f64 = self.λ[nb_λ - 1] - λ_1;
        let dλ  : f64 = self.λ[1] - λ_1;

        self.κ_c = Array1::<f64>::zeros(nb_λ);
        self.ϵ_c = Array1::<f64>::zeros(nb_λ);

        //let mut f : Array1<f64> = Array1::<f64>::zeros(nb_λ);
        let mut prof : VoigtProfile = VoigtProfile::new(&self.par, T, N);

        for iline in 0..nb_lines {
            let λ_ul = self.lines[iline].λ_ul;

            if λ_ul >= λ_1 && λ_ul <= λ_2 {
                let iλ0 = ((λ_ul - λ_1) / Δλ * ((nb_λ-1) as f64)) as i64;

                prof.set_Δλ(self.lines[iline].ΔλL);
                prof.set_mass(self.lines[iline].isotope_m);

                let diλ = std::cmp::max(2, ( (self.lines[iline].ΔλL + self.lines[iline].ΔλG) * self.par.Δλ_factor / dλ) as i64);
                let iλm = std::cmp::max(0, iλ0 - diλ) as usize;
                let iλp = std::cmp::min(iλ0 + diλ + 1, (nb_λ - 1) as i64) as usize;

                let κ = self.lines[iline].κ;
                let ϵ = self.lines[iline].ϵ;

                // paralellized
                //let t1 = Instant::now();
                //for iλ in im..ip {
                //    f[iλ] = prof.f(self.λ[iλ], λ_ul);
                //}
                //let t2 = Instant::now();

                //Zip::from(&mut self.κ_c.slice_mut(s![im..ip]))
                //    .and(&mut self.ϵ_c.slice_mut(s![im..ip]))
                //    .and(&f.slice(s![im..ip]))
                //    .par_apply(|a1, a2, &ff| {
                //        *a1 += κ * ff;
                //        *a2 += ϵ * ff;
                //    }
                //);
                //let t3 = Instant::now();

                for iλ in iλm..iλp {
                    let ff = prof.f(self.λ[iλ], λ_ul);
                    self.κ_c[iλ] =  self.κ_c[iλ] + κ * ff;
                    self.ϵ_c[iλ] =  self.ϵ_c[iλ] + ϵ * ff;
                }
            }
        }
    }

    /// integrate along the path
    pub fn integrate_along_path(&mut self, iN: usize, iθ: usize, NCO2: &f64, θ: &f64) -> ResultData {
        // z
        let nb_zsteps = self.z.len();
        let zp = self.p_vs_h.xy[[self.p_vs_h.xy.nrows()-1, 0]];
        let mut z_max  = self.t_vs_h.xy[[self.t_vs_h.xy.nrows()-1, 0]];
        if z_max > zp {
            z_max = zp;
        }
        if z_max > self.z[nb_zsteps-1] {
            z_max = self.z[nb_zsteps-1]
        }
        let mut znext = 10.0;
        let zoutmax = 1.2e4;
        let Dz = 2.0e3;

        let κds_limit = 0.01;

        // number of wavelengths and wavelength intervall
        let nb_λ = self.λ.len();
        let dλ = self.λ[1] - self.λ[0];

        // initial intensity
        let mut I_λ = self.planck_λ(self.par.T_surface, &self.λ);
        I_λ = I_λ * (1.0 - self.par.albedo);

        // pessure, temperature at z = 0
        let p = self.p_vs_h.y_at(0.0);
        let T = self.t_vs_h.y_at(0.0);
        let N = p / (self.par.kB * T);

        // results array
        let mut results = ResultData::new();
        let I_mean = I_λ.iter().sum::<f64>()*dλ;
        results.add(0.0, *NCO2, *θ, I_mean, T, N, 0.0, 0.0, 0.0, 0.0);

        // save results of zeroth step
        let out = format!("iz = {:3}, z = {:12.5e},  NCO2 = {:12.5e}, θ = {:12.5e}, I = {:12.5e}, T = {:12.5e}, N = {:12.5e}, ϵ = {:12.5e}, κ = {:12.5e}, ΔλL = {:12.5e}, ΔλD = {:12.5e}\n",
                                0, 0.0, NCO2*1.0e6, θ*180.0/PI, I_mean, T, N, 0.0, 0.0, 0.0, 0.0);
        match self.logfile.write_all(out.as_bytes()) {
            Err(why) => { panic!("couldn't write to logfile: {}", why.to_string()) },
            Ok(_) => 0
        };
        self.logfile.flush();
        print!("{}", out);

        // total emission
        let mut tot_e = 0.0;
        //
        let mut Tmin = self.par.T_surface;
        let mut Nmin = 1.0e30;

        for istep in 0..nb_zsteps {
            let mut vκ : Vec<f64> = Vec::<f64>::new();
            let mut vϵ : Vec<f64> = Vec::<f64>::new();

            // pressure, temperature and density at height = z
            let p = self.p_vs_h.y_at(self.z[istep]);
            let mut T = self.t_vs_h.y_at(self.z[istep]);
            let mut N = p / (self.par.kB * T);

            if self.par.T_of_h {
                T = self.par.T_surface;
            }
            if self.par.N_of_h {
                T = self.par.T_surface;
                N = p / (self.par.kB * T);
            }
            if T < Tmin {
                Tmin = T;
            }
            if N < Nmin {
                Nmin = N;
            }

            let t1 = Instant::now();
            // compute the line coefficients
            let (ΔλL_mean, ΔλD_mean) = self.compute_emission_and_absorption(T, N, p, *NCO2);
            let t2 = Instant::now();
            // compute ϵ(λ) = sum_i \int ϵ_i(λ_i) f(λ-λ_i) dλ
            // compute κ(λ) = sum_i \int κ_i(λ_i) f(λ-λ_i) dλ
            self.convolution(T, N);
            let t3 = Instant::now();

            // add background ?
            if self.par.background > 1.0e-10 {
                let iw = ((ΔλL_mean + ΔλD_mean) / dλ * self.par.Δλ_factor) as usize;
                // compute moving average
                let ma_κ = moving_average(&self.κ_c, iw*2);
                let ma_ϵ = moving_average(&self.ϵ_c, iw*2);
                // add background
                for i in 0..nb_λ {
                    self.κ_c[i] = self.κ_c[i] + ma_κ[i] * self.par.background;
                    self.ϵ_c[i] = self.ϵ_c[i] + ma_ϵ[i] * self.par.background;
                }
                // moving average once
                if iN == 0 && iθ == 0 {
                    self.save_intensity_as_npy(&self.λ, &ma_κ, [self.par.out_dir.as_str(), "moving_average_kappa"].join("/").as_str());
                }
            }

            // step size ds = z/cos(θ)
            let mut ds : f64 = 1.0;
            if istep < nb_zsteps - 1 {
                ds = (self.z[istep+1] - self.z[istep]) / θ.cos();
            }
            else {
                ds = (self.z[istep] - self.z[istep-1]) / θ.cos();
            }

            // integrate from s => s+ds at all wavelengths
            for i in 0..nb_λ {
                let κ = self.κ_c[i];
                let ϵ = self.ϵ_c[i];
                tot_e = tot_e + ϵ * ds;
                let exp_κ = (-κ * ds).exp();

                if self.par.with_emission {
                    let mut eps : f64 = 0.0;
                    if κ.abs() * ds < κds_limit {
                        eps = ϵ * ds;
                    }
                    else {
                        eps = ϵ / κ * (1.0 - exp_κ);
                    }
                    I_λ[i] = I_λ[i] * exp_κ + eps;
                }
                else {
                    I_λ[i] = I_λ[i] * exp_κ
                }
                vϵ.push(ϵ);
                vκ.push(I_λ[i] * κ);
            } // loop over λ

            // save intensity and spectrum
            if self.z_iout[istep] == 1 {
                znext = self.z[istep] + Dz;

                let Iname = format!("intensity_{}_{}_{:5.3e}", iN, iθ, self.z[istep]);
                self.save_intensity_as_npy(&self.λ, &I_λ, [self.par.out_dir.as_str(), &Iname].join("/").as_str());

                let cname = format!("spectrum_{}_{}_{:5.3e}", iN, iθ, self.z[istep]);
                let fname = [self.par.out_dir.as_str(), &cname].join("/");
                self.save_spectrum_as_npy(&fname);
            }

            // integrate over wavelength
            let I_mean = I_λ.iter().sum::<f64>()*dλ;
            let ϵ_mean = vϵ.into_iter().sum::<f64>()*dλ;
            let κ_mean = vκ.into_iter().sum::<f64>()*dλ;

            // set results vector
            results.add(self.z[istep], *NCO2, *θ, I_mean, T, N, ϵ_mean, κ_mean, ΔλL_mean, ΔλD_mean);

            // write results to log file
            let out = format!("iz = {:3}, z = {:12.5e},  NCO2 = {:12.5e}, θ = {:12.5e}, I = {:12.5e}, T = {:12.5e}, N = {:12.5e}, ϵ = {:12.5e}, κ = {:12.5e}, ΔλL = {:12.5e}, ΔλD = {:12.5e}\n",
                                    istep, self.z[istep], NCO2*1.0e6, θ*180.0/PI, I_mean, T, N, ϵ_mean, κ_mean, ΔλL_mean, ΔλD_mean);
            match self.logfile.write_all(out.as_bytes()) {
                Err(why) => { panic!("couldn't write to logfile: {}", why.to_string()) },
                Ok(_) => 0
            };
            self.logfile.flush();
            print!("{}", out);
            //self.resultdata.add(self.z[istep], NCO2*1.0e6, θ*180.0/PI, I_mean, T, N, ϵ_mean, κ_mean, ΔλL_mean, ΔλD_mean);
        }  // lop over z istep

        results
    }

    /// Compute absorption
    pub fn integrate(&mut self) {

        // compute Planck intensity
        let λ1 = 1.0e-6;
        let λ2 = 1.0e-4;
        let nλ : usize = 1000;
        let T = self.par.T_surface;
        self.compute_and_save_planck(λ1, λ2, nλ, T, [self.par.out_dir.as_str(), "planck_intensity"].join("/").as_str());

        // vector of angles
        let vθ = &self.theta_deg * PI/180.0;
        // vector of CO2 concentrations
        let vNCO2 = self.NCO2.clone();

        //logfile header
        let out = format!("{:3} {:3} {:3} {:3}\n", self.NCO2.len(), vθ.len(), self.z.len(), 10);
        match self.logfile.write_all(out.as_bytes()) {
            Err(why) => { panic!("couldn't write to logfile: {}", why.to_string()) },
            Ok(_) => 0
        };
        print!("{}", out);

        // loop over CO2 concentrations
        for (iN, NCO2) in vNCO2.iter().enumerate()  {
            // loop over angles
            for (iθ, θ) in vθ.iter().enumerate()  {

                // intermediate log file header line
                let out = format!("# {:3} {:3} {:12.5e} {:12.5e}\n", iN, iθ, (*NCO2)*1.0e6, θ / PI * 180.0);
                match self.logfile.write_all(out.as_bytes()) {
                    Err(why) => { panic!("couldn't write to logfile: {}", why.to_string()) },
                    Ok(_) => 0
                };
                print!("{}", out);

                // integrate along path
                let mut results = self.integrate_along_path(iN, iθ, NCO2, θ);

                // save results
                let fname = format!("result-{}_{}.npy", iN, iθ);
                results.write_data(&self.par.out_dir.as_str(), &fname);
            }
        }
        //self.resultdata.write_data(self.par.out_dir.as_str(), "data.npy");
    }
}