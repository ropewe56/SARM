// https://pyo3.rs/v0.7.0/function.html?highlight=pyfn#python-functions
// https://github.com/Dushistov/rust_swig
// https://github.com/rust-lang/rust-bindgen
// https://lib.rs/development-tools/ffi

#![feature(non_ascii_idents)]
#![allow(non_snake_case)]
#![allow(unused)]

extern crate libsarm;
use libsarm::spectrum;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyfunction]
fn call_spectrum(json_file: &str) -> PyResult<String> {
    let mut spec = libsarm::spectrum::spectrum::Spectrum::new(&json_file);
    spec.compute_all();
    Ok("end".to_string())
}

#[pyfunction]
fn other(json_file: &str) -> PyResult<String> {
    let mut spec = libsarm::spectrum::spectrum::Spectrum::new(&json_file);
    spec.compute_all();
    Ok("end".to_string())
}

#[pyfunction]
fn sum_as_string_py(_py: Python, a:i64, b:i64) -> PyResult<String> {
    let out = sum_as_string(a, b);
    Ok(out)
}

// This module is a python module implemented in Rust.
#[pymodule]
fn pysarm(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(call_spectrum))?;
    m.add_wrapped(wrap_pyfunction!(other))?;
    m.add_wrapped(wrap_pyfunction!(sum_as_string_py))?;

    Ok(())
}

//// add bindings to the generated python module
//// N.B: names: "librust2py" must be the name of the `.so` or `.pyd` file
///// This module is implemented in Rust.
//#[pymodule]
//fn rust2py(py: Python, m: &PyModule) -> PyResult<()> {
//
//    // PyO3 aware function. All of our python interfaces could be declared in a separate module.
//    // Note that the `#[pyfn()]` annotation automatically converts the arguments from
//    // Python objects to Rust values; and the Rust return value back into a Python object.
//    #[pyfn(m, "sum_as_string")]
//    fn sum_as_string_py(_py: Python, a:i64, b:i64) -> PyResult<String> {
//       let out = sum_as_string(a, b);
//       Ok(out)
//    }
//
//    Ok(())
//}

// logic implemented as a normal rust function
fn sum_as_string(a:i64, b:i64) -> String {
    format!("{}", a + b).to_string()
}
