use ndarray::prelude::*;
use std::f64::consts::PI;

use crate::spectrum::parameter::*;

pub struct GaussProfile {
    pub mass : f64,
    pub kB   : f64,
    pub c    : f64,
    pub π    : f64,
    pub ln2  : f64,
    pub temp : f64,
    pub Δλ0  : f64,
    pub b    : f64,
}

impl GaussProfile {
    pub fn new(par: &Parameter, temp: f64) -> GaussProfile {
        let ln2 = (2.0 as f64).ln();
        let Δλ0 = (2.0 * par.kB * temp / par.mCO2).sqrt() / par.c0;
        let b   = (1.0 / PI).sqrt();
        GaussProfile{mass: par.mCO2, kB: par.kB, c: par.c0, π : PI, ln2: ln2, temp: temp,  Δλ0: Δλ0, b: b}
    }

    pub fn set_T(&self, temp: f64) {
        let Δλ0 = (2.0 * self.kB * temp / self.mass).sqrt() / self.c;
        let b   = (1.0 / self.π).sqrt();
    }

    pub fn set_mass(&mut self, mass: f64) {
        self.mass = mass;
    }

    #[inline]
    pub fn f(&self, λ: f64, λ0: f64) -> f64 {
        let Δλ = self.Δλ0 * λ0;
        //println!("{:12.6e}  {:12.6e}", self.Δλ0, Δλ);
        let c  = self.ln2 / (Δλ*Δλ);
        self.b * c.sqrt() * (- c * (λ - λ0)*(λ - λ0)).exp()
    }

    pub fn fa(&self, λ: &Array1<f64>, λ0: f64) -> Array1<f64> {
        let Δλ = self.Δλ0 * λ0;
        let c  = self.ln2 / (Δλ*Δλ);
        let mut v : Array1<f64> = (- c *(λ - λ0)*(λ - λ0));
        v.mapv_inplace(|x| x.exp());
        self.b * v
    }
}

pub struct LorentzProfile {
    Δλ : f64,
    π  : f64,
}

impl LorentzProfile {
    pub fn new(par: &Parameter) -> LorentzProfile {
        LorentzProfile{Δλ: 1.0e-9, π: PI}
    }

    #[inline]
    pub fn f(&self, λ: f64, λ0: f64) -> f64 {
        let dλ = (λ - λ0)/self.Δλ;
        1.0 / (self.π * self.Δλ * (dλ*dλ + 1.0))
    }

    #[inline]
    pub fn fa(&self, λ: &Array1<f64>, λ0: f64) -> Array1<f64> {
        let mut dλ = (λ - λ0)/self.Δλ;
        dλ.mapv_inplace(|x| x*x);
        1.0 / (self.π * self.Δλ * (dλ + 1.0))
    }

    pub fn set_Δλ(&mut self, Δλ: f64) {
        self.Δλ = Δλ;
    }
}

pub struct VoigtProfile {
    Gauss   : GaussProfile,
    Lorentz : LorentzProfile,
}

impl VoigtProfile {
    pub fn new(par: &Parameter, temp: f64, density: f64) -> VoigtProfile {
        let gauss   = GaussProfile::new(par, temp);
        let lorentz = LorentzProfile::new(par);

        VoigtProfile{Gauss: gauss, Lorentz: lorentz}
    }

    pub fn set_Δλ(&mut self, ΔλL: f64) {
        self.Lorentz.set_Δλ(ΔλL);
    }

    pub fn set_mass(&mut self, mass: f64) {
        self.Gauss.set_mass(mass);
    }

    #[inline]
    pub fn f(&self, λ: f64, λ0: f64) -> f64 {
        let mut v = (self.Lorentz.Δλ / (self.Gauss.Δλ0*λ0));
        v = 1.36606 * v - 0.47719 *v*v + 0.11116 * v*v*v;
        if v > 1.0 {
            let nix = self.Gauss.f(λ, λ0);
            self.Lorentz.f(λ, λ0)
        }
        else {
            v * self.Lorentz.f(λ, λ0) + (1.0 - v) * self.Gauss.f(λ, λ0)
        }
    }

    #[inline]
    pub fn fa(&self, λ: &Array1<f64>, λ0: f64) -> Array1<f64> {
        let mut v = (self.Lorentz.Δλ / (self.Gauss.Δλ0*λ0));
        v = 1.36606 * v - 0.47719 *v*v + 0.11116 * v*v*v;
        if v > 1.0 {
            self.Lorentz.fa(λ, λ0)
        }
        else {
            v * self.Lorentz.fa(λ, λ0) + (1.0 - v) * self.Gauss.fa(λ, λ0)
        }
    }


//    """
//    Return the Voigt line shape at x with Lorentzian component HWHM gamma
//    and Gaussian component HWHM alpha.
//
//    """
//    sigma = alpha / np.sqrt(2 * np.log(2))
//
//    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma /np.sqrt(2*np.pi)

// wofz Faddeeva function
//
// Returns the value of the Faddeeva function for complex argument:
//
// exp(-z**2) * erfc(-i*z)


}
