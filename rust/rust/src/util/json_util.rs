use std::io;
use std::io::prelude::*;
use std::fs::File;
use serde_json::{Value};

pub fn read_json(jsonfname: &str) -> Result<Value, io::Error> {
    let mut f = File::open(jsonfname)?;
    let mut buffer = String::new();
    f.read_to_string(&mut buffer)?;
    let v: Value = serde_json::from_str(&buffer)?;
    Ok(v)
}
