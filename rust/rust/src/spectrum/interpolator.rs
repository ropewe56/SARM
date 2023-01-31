use ndarray::prelude::*;
use crate::spectrum::ioutil::{save_npy_2, save_npy_1, load_npy_2, load_npy_1, open_file};

pub struct Interpolator {
    pub xy  : Array2<f64>,
    pub x1  : f64,
    pub ddx : f64,
}

impl Interpolator {

    pub fn new(fname: &str) -> Interpolator {
        let xy = load_npy_2(&fname);
        println!("{} {} {}", xy.nrows(), xy.ncols(), fname);
        let n : usize = xy.nrows();
        let ddx = 1.0 / (xy[[1,0]] - xy[[0,0]]);
        let x1  = xy[[0,0]];

        Interpolator{xy: xy, x1: x1, ddx: ddx}
    }

    pub fn new_array(xy0: &Array2<f64>) -> Interpolator {
        let xy = xy0.clone();
        let n : usize = xy.nrows();
        let ddx = 1.0 / (xy[[1,0]] - xy[[0,0]]);
        let x1  = xy[[0,0]];

        Interpolator{xy: xy, x1: x1, ddx: ddx}
    }

    pub fn y_at(&self, x: f64) -> f64 {
        let ix : usize = std::cmp::min( ((x - self.x1) * self.ddx) as usize, self.xy.nrows()-2 );
        let u = (x - self.xy[[ix,0]]) * self.ddx;
        (1.0 - u) * self.xy[[ix,1]] + u * self.xy[[ix+1,1]]
    }

    pub fn yy_at(&self, x: f64, id: usize) -> f64 {
        let ix : usize = std::cmp::min( ((x - self.x1) * self.ddx) as usize, self.xy.nrows()-2 );
        let u = (x - self.xy[[ix,0]]) * self.ddx;
        if self.xy.ncols() <= id+1 {
            println!("{} {}", self.xy.ncols(), id)
        }
        (1.0 - u) * self.xy[[ix,id+1]] + u * self.xy[[ix+1,id+1]]
    }
}
