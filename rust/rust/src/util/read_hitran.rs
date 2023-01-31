use std::io;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufReader;
//use std::str::FromStr;
//use nalgebra::{Vector3};
//type Vector3d = Vector3<f64>;

use ndarray::prelude::*;
//use std::ops::Mul;
//use std::ops::Div;
//use std::ops::Add;
//use std::ops::Sub;
//use ndarray::Data;
//use ndarray::Dimension;
//use ndarray::Shape;

//type Array1d = Array1<f64>;
//type Array2d = Array2<f64>;

use crate::util::hdf5_util;
use crate::util::parameter::Parameter;

pub struct Hitran {
    pub width: Vec<usize>,
    pub pos: Vec<usize>,
    pub h: f64,
    pub c: f64,
}

impl Hitran {

    pub fn new(par: &Parameter) -> Hitran {
        let width1 : Vec<usize> = vec![  2,   1,  12,  10,  10,   5,   5,  10,   4,   8,  15,  15,  15,  15];
        let width2 : Vec<usize> = vec![ 1,  2,  1,  7,  7];
        let w1 : usize = width1.iter().sum();
        let w2 : usize = width2.iter().sum();
        let n = 160 as usize - w1 - w2;
        let width = vec![  2,   1,  12,  10,  10,   5,   5,  10,   4,   8,  15,  15,  15,  15, n,  1,  2,  2,  7,  7];

        let mut pos = Vec::new();
        let mut sum_of_elems: usize = 0;
        for n in 0..width.len() {
            pos.push(sum_of_elems);
            sum_of_elems += width[n];
        }
        let j: usize = pos.len();
        pos[j-1] = pos[j-1] - 1;

        Hitran{width: width, pos: pos, h: par.h, c: par.c}
    }

    fn get_float(&self, s: &str, i: usize) -> f64 {
        let s_ = (&s[self.pos[i]..self.pos[i] + self.width[i]]).trim();
        //println!("{} |{}| {} {}", i, s_, self.pos[i], self.width[i]);
        let f = s_.parse::<f64>().unwrap();
        f
    }

    fn get_int(&self, s: &str, i: usize) -> i32 {
        let s_ = (&s[self.pos[i]..self.pos[i] + self.width[i]]).trim();
        //println!("{} |{}| {} {}", i, s_, self.pos[i], self.width[i]);
        let f = s_.parse::<i32>().unwrap();
        f
    }

    //fn get_str(&self, s: &str, i: usize) -> String {
    //    let s_ = (&s[self.pos[i]..self.pos[i] + self.width[i]]).trim();
    //    String::new(s_)
    //}


    pub fn extract_line_par(&self, line: &String) -> Vec<f64> {
        println!("{}", line.len());
        let mid    = self.get_int(&line,  0);
        let iid    = self.get_int(&line,  1);
        let ν      = self.get_float(&line,  2);
        let S      = self.get_float(&line,  3);
        let A      = self.get_float(&line,  4);
        let γ_a    = self.get_float(&line,  5);
        let γ_s    = self.get_float(&line,  6);
        let ν_l    = self.get_float(&line,  7);
        let n_a    = self.get_float(&line,  8);
        let δ_a    = self.get_float(&line,  9);
        //let guq    = self.get_str(&line, 10);
        //let glq    = self.get_str(&line, 11);
        //let luq    = self.get_str(&line, 12);
        //let llq    = self.get_str(&line, 13);
        let error  = self.get_int(&line, 15);
        let refer  = self.get_int(&line, 16);
        //let lmflag = self.get_str(&line, 17);
        let g_u    = self.get_float(&line, 18);
        let g_l    = self.get_float(&line, 19);


        let ν      = ν   / 1.0e-2;  // 1/m
        let ν_l    = ν_l / 1.0e-2;  // 1/m
        let γ_a    = γ_a / 1.0e-2;  // 1/m
        let γ_s    = γ_s / 1.0e-2;  // 1/m

        let λ      = 1.0 / ν;
        let ΔE_ul  = self.h * self.c / λ;
        let E_l    = self.h * self.c * ν_l;
        let E_u    = E_l + ΔE_ul;
        let Δλa    = γ_a * λ * λ; // hwhm
        let Δλs    = γ_s * λ * λ; // hwhm

        let mut c = vec![0.0; 13];

        // S = 1/N * (N_l * B_lu - N_u * B_ul) * h * ν / c
        // A_ul = 8 * π * h * ν**3 * B_ul
        // g_l * B_lu = g_u * B_ul
        // S_ν^N(T) = 1/Q_{tot}(T) * ( exp(-E_l/(kB*T)) * gl * B_lu - exp(-E_u/(kB*T)) * gu * B_ul) * h * ν_0 / c
        // S_ν^N(T) = g_1/Q_{tot}(T) * A_ul / (8 π c \nu_0**2) exp(-E_l/(kB*T)) * (1.0 - epx(-h ν / (kB*T)))

        c[ 0] = λ;    // m
        c[ 1] = E_l;  // eV
        c[ 2] = E_u;  // eV
        c[ 3] = S;    // 1/s
        c[ 4] = A;    // 1/s
        c[ 5] = γ_a * 1.0e-5;  // 1/m/pascal
        c[ 6] = γ_s * 1.0e-5;  // 1/m/pascal
        c[ 7] = n_a;  //
        c[ 8] = δ_a / 1.0e-2 * 1.0e-5;  // δ_a : [1/cm 1/atm] => [1/1e-2m 1/10e5p] => 1/m * 1/pascal

        return c;
    }

    pub fn read_par_par(&self, hitran_name: &str) -> Result<Array2<f64>, io::Error> {
        let mut file = File::open(hitran_name)?;
        let mut vv = Vec::new();
        for line in BufReader::new(file).lines() {
            let vl: Vec<f64> = self.extract_line_par(&(line?));
            vv.push(vl);
        }
        Ok(crate::util::hdf5_util::vecvec_to_array2(&vv))
    }

    pub fn extract_line(&self, par: &Parameter, v: Vec<&str>) -> Option<Vec<f64>> {

        let mid = v[ 0].parse::<f64>().unwrap();
        let iid = v[ 1].parse::<f64>().unwrap();
        let ν   = v[ 2].parse::<f64>().unwrap();
        let ν   = ν   / 1.0e-2;  // 1/m
        let λ   = 1.0 / ν;

        if par.λmin <= λ && λ <= par.λmax {

            let S   = v[ 3].parse::<f64>().unwrap();
            let A   = v[ 4].parse::<f64>().unwrap();
            let γ_a = v[ 5].parse::<f64>().unwrap();
            let γ_s = v[ 6].parse::<f64>().unwrap();
            let ν_l = v[ 7].parse::<f64>().unwrap();
            let n_a = v[ 8].parse::<f64>().unwrap();
            let δ_a = v[ 9].parse::<f64>().unwrap();
            let g_u = v[10].parse::<f64>().unwrap();
            let g_l = v[11].parse::<f64>().unwrap();

            let ν_l    = ν_l / 1.0e-2;  // 1/m
            let γ_a    = γ_a / 1.0e-2 * 1.0e-5;  // 1 / (m pascal)
            let γ_s    = γ_s / 1.0e-2 * 1.0e-5;  // 1 / (m pascal)
            let δ_a    = δ_a / 1.0e-2 * 1.0e-5;  //
            let ΔE_ul  = self.h * self.c / λ;
            let E_l    = self.h * self.c * ν_l;
            let E_u    = E_l + ΔE_ul;

            let mut c = vec![0.0; 13];

            // S = 1/N * (N_l * B_lu - N_u * B_ul) * h * ν / c
            // A_ul = 8 * π * h * ν**3 * B_ul
            // g_l * B_lu = g_u * B_ul
            // S_ν^N(T) = 1/Q_{tot}(T) * ( exp(-E_l/(kB*T)) * gl * B_lu - exp(-E_u/(kB*T)) * gu * B_ul) * h * ν_0 / c
            // S_ν^N(T) = g_1/Q_{tot}(T) * A_ul / (8 π c \nu_0**2) exp(-E_l/(kB*T)) * (1.0 - epx(-h ν / (kB*T)))

            c[ 0] = λ;    // m
            c[ 1] = E_l;  // eV
            c[ 2] = E_u;  // eV
            c[ 3] = S;    // 1/s
            c[ 4] = A;    // 1/s
            c[ 5] = γ_a;  // 1 / (m pascal)
            c[ 6] = γ_s;  // 1 / (m pascal)
            c[ 7] = n_a;  //
            c[ 8] = δ_a;  //  1 / (m pascal)
            c[ 9] = g_u;  //
            c[10] = g_l;  //

            Some(c)
        }
        else {
            None
        }
    }

    pub fn read_par(&self, par: &Parameter, hitran_name: &str) -> Result<Array2<f64>, io::Error> {
        let mut file = File::open(hitran_name)?;
        let mut vv = Vec::new();
        let mut i = 0;
        let mut smax = 0.0;
        let num_header_lines = 1;
        for line in BufReader::new(file).lines() {
            i = i + 1;
            if i > num_header_lines {
                let ll = line?;
                let mut split = ll.split(",");
                let vec: Vec<&str> = split.collect();

                let vl = self.extract_line(par, vec);
                match vl {
                    Some(x) => {vv.push(x); smax = if smax < x[3] {x[3]} else {smax}}, 
                    None => (),
                }
            }
        }
        println!("smax{}", smax);
        Ok(crate::util::hdf5_util::vecvec_to_array2(&vv))
    }
}
