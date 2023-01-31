use ndarray::prelude::*;
//use rayon::prelude::*;
//use ndarray_parallel::prelude::*;
//use ndarray::Zip;
use ndarray::{Data, RemoveAxis, Zip};

use std::cmp;
use std::time::{Duration, Instant};
use std::fmt::format;

use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
//use std::path::Path;
//use std::path::PathBuf;

use std::collections::HashMap;

use crate::spectrum::spectrum::Spectrum;

use std::f64::consts::PI;

impl Spectrum {
    /// integrate F (Reinhart) along path 
    pub fn integrate_along_path_F(&mut self, iN: usize, iθ: usize, NCO2: &f64, θ: &f64) {

        let nz = self.z.len();
        let zp = self.p_vs_h.xy[[self.p_vs_h.xy.nrows()-1, 0]];
        let mut z_max  = self.t_vs_h.xy[[self.t_vs_h.xy.nrows()-1, 0]];
        if z_max > zp {
            z_max = zp;
        }
        if z_max > self.z[nz-1] {
            z_max = self.z[nz-1]
        }

        let mut Tmin = self.par.T_surface;
        let mut Nmin = 1.0e30;

        let nb_steps = nz;

        let p = self.p_vs_h.y_at(0.0);
        let T = self.t_vs_h.y_at(0.0);
        let N = p / (self.par.kB * T);

        let nb_lines = self.lines.len();

        let mut znext = 10.0;
        let zoutmax = 1.2e4;
        let Dz = 2.0e3;

        for iline in 0..nb_lines {
            self.lines[iline].κ_sum = 0.0;
        }

        for istep in 0..nb_steps {

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

            let (ΔλL_mean, ΔλD_mean) = self.compute_emission_and_absorption(T, N, p, *NCO2);

            let mut ds : f64 = 1.0;
            if istep < nb_steps-1 {
                ds = (self.z[istep+1] - self.z[istep]) / θ.cos();
            }
            else {
                ds = (self.z[istep] - self.z[istep-1]) / θ.cos();
            }

            for iline in 0..nb_lines {
                self.lines[iline].κ_sum += self.lines[iline].κ * ds;
            }
        }  // istep
    }

    /// Compute F (Reinhart) 
    pub fn compute_F(&mut self) {

        let nb_lines  = self.lines.len();
        let nb_λ      = self.λ.len();
        let λ_1 : f64 = self.λ[0];
        let λ_2 : f64 = self.λ[nb_λ - 1];
        let Δλ  : f64 = self.λ[nb_λ - 1] - λ_1;
        let dλ  : f64 = self.λ[1] - λ_1;
        let I_λ       = self.planck_λ(self.par.T_surface, &self.λ);
        let I_mean    = I_λ.iter().sum::<f64>()*dλ;


        // vector of angles
        let vθ = &self.theta_deg * PI/180.0;
        // vector of CO2 concentraions
        let vNCO2 = self.NCO2.clone();

        // Stefan-Boltzmann
        let S0 = self.par.σ * f64::powf(self.par.T_surface, 4.0);

        // not used
        for iline in 0..nb_lines {
            self.lines[iline].f_u = self.planck(self.par.T_surface, self.lines[iline].λ_ul);
        }

        // loop over CO2 cocentrations
        for (iN, NCO2) in vNCO2.iter().enumerate()  {
            // loop over angles
            for (iθ, θ) in vθ.iter().enumerate()  {

                self.integrate_along_path_F(iN, iθ, NCO2, θ);

                let mut e_sum1 = 0.0;
                let mut e_sum2 = 0.0;
                for iline in 0..nb_lines {
                    e_sum1 += self.lines[iline].κ_sum;
                    e_sum2 += 1.0 - (-self.lines[iline].κ_sum / self.lines[iline].ΔλL).exp();
                }

                let e_sum_exp = (- e_sum1 / Δλ).exp();

                let out = format!("# {:3} {:3}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}  {:12.5e}\n",
                    iN, iθ, (*NCO2)*1.0e6, θ / PI * 180.0,
                    I_mean, I_mean*e_sum_exp, e_sum1, e_sum2);

                match self.logfile_F.write_all(out.as_bytes()) {
                    Err(why) => { panic!("couldn't write to logfile_F: {}", why.to_string()) },
                    Ok(_) => 0
                };
                print!("{}", out);
                //save_npy_1("e_sum.npy", &e_sum);
            }
        }
    }
}