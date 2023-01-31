// https://pyo3.rs/v0.7.0/function.html?highlight=pyfn#python-functions
// https://github.com/Dushistov/rust_swig
// https://github.com/rust-lang/rust-bindgen
// https://lib.rs/development-tools/ffi

#![allow(non_snake_case)]
#![allow(unused)]

extern crate sarm;
use sarm::spectrum;

use numpy::pyo3::Python;
use numpy::ndarray::array;
use numpy::{ToPyArray, PyArray};
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyfunction]
fn call_spectrum(json_file: &str) -> PyResult<String> {
    let mut spec = sarm::spectrum::spectrum::Spectrum::new(&json_file);
    spec.integrate();
    Ok("end".to_string())
}


// This module is a python module implemented in Rust.
#[pymodule]
fn pysarm(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(call_spectrum))?;
    Ok(())
}
