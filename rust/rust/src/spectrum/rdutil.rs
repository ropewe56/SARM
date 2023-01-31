use ndarray::prelude::*;

#[inline]
pub fn float_max(a: &Vec<f64>) -> f64 {
    let mut b = -1.66e66;
    for i in 0..a.len() {
        if a[i] >= b {
            b = a[i];
        }
    }
    b
}

#[inline]
pub fn float_min(a: &Vec<f64>) -> f64 {
    let mut b = 1.66e66;
    for i in 0..a.len() {
        if a[i] <= b {
            b = a[i];
        }
    }
    b
}

#[inline]
pub fn float_mean(a: &Vec<f64>) -> f64 {
    let asum : f64 = a.iter().sum();
    asum / (a.len() as f64)
}

pub fn moving_average(y: &Array1<f64>, iw: usize) -> Array1<f64> {
    let n = y.len();
    let mut ma : Array1<f64> = Array1::<f64>::zeros(n);
    let mut mi : Array1<f64> = Array1::<f64>::zeros(n);
    ma[0] = y.slice(s![0..iw]).sum();
    mi[0] = iw as f64;
    let mut ym = 0.0;
    let mut mam = 0.0;
    let mut mim = 0.0;

    for i in 1..n {

        let im = (i as i64) - (iw as i64);
        let ip = (i as i64) + (iw as i64);

        let mut dma1 = 0.0;
        let mut dmi1 = 0.0;
        let mut dma2 = 0.0;
        let mut dmi2 = 0.0;

        if im >= 0 {
            dma1 = - y[im as usize];
            dmi1 = - 1.0;
        }
        if ip <= (n-1) as i64 {
            dma2 = dma2 + y[ip as usize];
            dmi2 = dmi2 + 1.0;
        }
        ma[i] = ma[i-1] + dma1 + dma2;
        mi[i] = mi[i-1] + dmi1 + dmi2;

        ym += y[i];
    }

    for i in 0..n {
        if mi[i] > 0.0 {
            ma[i] = ma[i] / mi[i];
        }
        mam += ma[i];
        mim += mi[i];
    }

    //save_npy_1(["/home/wester/Projects/simpleCO2/radoutput/e/00_0_1_1.0_0_0.01_1e-11", "ma_k0.npy"].join("/").as_str(), &ma);
    //save_npy_1(["/home/wester/Projects/simpleCO2/radoutput/e/00_0_1_1.0_0_0.01_1e-11", "mi_k0.npy"].join("/").as_str(), &mi);
    //println!("{:12.6e}  {:12.6e}  {:12.6e}  n = {} iw = {}", ym, mam, mim, n, iw);

    ma
}
