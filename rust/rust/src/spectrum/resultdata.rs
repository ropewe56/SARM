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
use std::path::Path;
use std::path::PathBuf;

use crate::spectrum::ioutil::{save_npy_2};

pub struct ResultData {
    pub data : Vec<Vec<f64>>,
}

impl ResultData {
    pub fn new() -> ResultData {
        let mut data : Vec<Vec<f64>> = Vec::<Vec<f64>>::new();
        ResultData{data : data}
    }

    pub fn add(&mut self, a1: f64, a2: f64, a3: f64, a4: f64, a5: f64, a6: f64, a7: f64, a8: f64, a9: f64, a10: f64) {
        let v = vec![a1, a2, a3, a4, a4, a6, a7, a8, a9, a10];
        self.data.push(v);
    }

    pub fn write_data(&mut self, dir: &str, fname: &str) {   
        let mut path = PathBuf::new();
        path.push(dir);
        path.push(fname);

        let n = self.data.len();
        let mut a = Array2::<f64>::zeros((n, 10));
        for i in 0..n {
            for j in 0..10 {
                a[[i,j]] = self.data[i][j];
            }
        }

        let pp : Option<&str> = path.to_str();
        match pp {
            Some(p) => { save_npy_2(p, &a) },
            None    => Ok(())
        };        
    }
}
