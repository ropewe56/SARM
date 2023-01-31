//use std::io;https://www.forrestthewoods.com/blog/how-to-debug-rust-with-visual-studio-code/

#![allow(confusable_idents)]
#![allow(uncommon_codepoints)]
#![allow(non_snake_case)]
#![allow(unused)]

//extern crate clap;
//use clap::{Arg, Command};
extern crate dirs;
extern crate physical_constants;
use std::path::Path;
use std::path::PathBuf;
use std::env;

extern crate sarm;
use sarm::spectrum::spectrum;

fn main() { //  -> hdf5::Result<()>)

    let args: Vec<String> = env::args().collect();
    //dbg!(&args);

    let json_file = args[1].as_str();

    let mut spec = sarm::spectrum::spectrum::Spectrum::new(&json_file);
    spec.integrate();
}