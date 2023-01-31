#![allow(uncommon_codepoints)]
#![allow(confusable_idents)]
#![allow(non_snake_case)]
#![allow(unused)]


#[macro_use]
extern crate ndarray;
extern crate ndarray_parallel;
extern crate ndarray_npy;

extern crate rayon;
use rayon::prelude::*;

use ndarray::prelude::*;
use ndarray_parallel::prelude::*;

//#[macro_use]
//use std::io::prelude::*;

extern crate clap;
extern crate dirs;
use std::path::Path;
use std::path::PathBuf;
use std::env;

pub mod spectrum;