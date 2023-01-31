use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::path::PathBuf;
use ndarray::prelude::*;
use ndarray_npy::{ReadNpyExt, WriteNpyError, ReadNpyError, WriteNpyExt};
use std::error::Error;

pub fn save_npy_2(fname: &str, arr: &Array2<f64>) -> Result<(), WriteNpyError> {
    let mut path = PathBuf::new();
    path.push(fname);
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };
    arr.write_npy(file)?;
    Ok(())
}

pub fn save_npy_1(fname: &str, arr: &Array1<f64>) -> Result<(), WriteNpyError> {
    let mut path = PathBuf::new();
    path.push(fname);
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };
    arr.write_npy(file)?;
    Ok(())
}

pub fn load_npy_2(fname: &str) -> Array2<f64> {
    let mut path = PathBuf::new();
    path.push(fname);
    let mut file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };
    let reader = File::open(path).unwrap();
    let arr = Array2::<f64>::read_npy(reader).unwrap();
    arr
}

pub fn load_npy_1(fname: &str) -> Array1<f64> {
    let mut path = PathBuf::new();
    path.push(fname);
    let mut file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };
    let reader = File::open(path).unwrap();
    let arr = Array1::<f64>::read_npy(reader).unwrap();
    arr
}

pub fn load_npy_1_i(fname: &str) -> Array1<i64> {
    let mut path = PathBuf::new();
    path.push(fname);
    let mut file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };
    let reader = File::open(path).unwrap();
    let arr = Array1::<i64>::read_npy(reader).unwrap();
    arr
}

pub fn open_file(dir: &str, fname: &str) -> File {
    let mut path = PathBuf::new();
    path.push(dir);
    path.push(fname);
    let mut logfile = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why.to_string()),
        Ok(logfile) => logfile,
    };
    logfile
}
