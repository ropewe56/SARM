use ndarray::prelude::*;

pub fn vecvec_to_array2(vv: &Vec<Vec<f64>>) -> Array2<f64> {
    let n1 = vv.len();
    let n2 = vv[0].len();

    let mut aa = Array2::<f64>::zeros((n1, n2));

    for i in 0..n1
    {
        for j in 0..n2
        {
            aa[[i, j]] = vv[i][j];
        }
    }
    aa
}

pub fn save_hdf5_1_V(hdf5_name: &str, aa: &Vec<f64>, grname: &str, dsname: &str) {
    let n1    = aa.len();

    let mut file : hdf5::File;
    let file_ok = hdf5::File::open(hdf5_name, "w");
    match file_ok {
        Ok(file_ok) => { file = file_ok; }
        Err(err) => { println!("error: {}. Could not open {}", err, hdf5_name); std::process::exit(1); }
    };


    let group = file.create_group(grname).unwrap();
    let data  = group.new_dataset::<f64>().create(dsname, n1).unwrap();
    data.write(aa).unwrap();
}

pub fn save_hdf5_1(hdf5_name: &str, aa: &Array1::<f64>, grname: &str, dsname: &str) {
    let n1    = aa.len();

    let mut file : hdf5::File;
    let file_ok = hdf5::File::open(hdf5_name, "w");
    match file_ok {
        Ok(file_ok) => { file = file_ok; }
        Err(err) => { println!("error: {}. Could not open {}", err, hdf5_name); std::process::exit(1); }
    };


    let group = file.create_group(grname).unwrap();
    let data  = group.new_dataset::<f64>().create(dsname, n1).unwrap();
    data.write(aa).unwrap();
}

pub fn load_hdf5_1(hdf5_name: &str, grname: &str, dsname: &str) -> Array1<f64> {
    let file = hdf5::File::open(hdf5_name, "r").unwrap();
    let dataset = file.dataset(&[grname, "/", dsname].concat()).unwrap();
    let data = dataset.read_1d().unwrap();
    data
}

pub fn save_hdf5_2(hdf5_name: &str, aa: &Array2::<f64>, grname: &str, dsname: &str) {
    let n1 = aa.nrows();
    let n2 = aa.ncols();
    let file = hdf5::File::open(hdf5_name, "w").unwrap();
    let group = file.create_group(grname).unwrap();
    let data = group.new_dataset::<f64>().create(dsname, (n1, n2)).unwrap();
    data.write(aa).unwrap();
}

pub fn save_hdf5_V2(hdf5_name: &str, aa: &Vec<&Array1::<f64>>, grname: &str, dsname: &str) {
    let n1 = aa.len();
    let n2 = aa[0].len();
    let mut A = Array2::<f64>::zeros((n1, n2));
    for i in 0..n2 {
        for j in 0..n1 {
            A[[i, j]] = aa[j][i];
        }
    }
    save_hdf5_2(hdf5_name, &A, grname, dsname);
}

pub fn load_hdf5_2(hdf5_name: &str, grname: &str, dsname: &str) -> Array2<f64> {
    let file = hdf5::File::open(hdf5_name, "r").unwrap();
    let dataset = file.dataset(&[grname, "/", dsname].concat()).unwrap();
    let data = dataset.read_2d().unwrap();
    data
}
