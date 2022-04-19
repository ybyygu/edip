// [[file:../edip.note::1141962a][1141962a]]
//#![deny(warnings)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
use gut::prelude::*;

type Array3 = [f64; 3];
// 1141962a ends here

// [[file:../edip.note::1ace2574][1ace2574]]
mod edip;
// 1ace2574 ends here

// [[file:../edip.note::f46cb3ad][f46cb3ad]]
// tau(Z) (Ismail & Kaxiras, 1993)
const u1: f64 = -0.165799;
const u2: f64 = 32.557;
const u3: f64 = 0.286198;
const u4: f64 = 0.66;

/// Parameters for EDIP
#[derive(Debug, Copy, Clone)]
pub struct EdipParams {
    A: f64,
    B: f64,
    rh: f64,
    a: f64,
    sig: f64,
    lam: f64,
    gam: f64,
    b: f64,
    c: f64,
    mu: f64,
    Qo: f64,
    bet: f64,
    alp: f64,
    /// cutoff for g(r)
    bg: f64,
    /// justo prefactor for bond order
    palp: f64,
    delta: f64,
    eta: f64,
}

impl EdipParams {
    /// EDIP-Si Parameters
    pub fn silicon() -> Self {
        let A = 5.6714030;
        let B = 2.0002804;
        let rh = 1.2085196;
        let a = 3.1213820;
        let sig = 0.5774108;
        let lam = 1.4533108;
        let gam = 1.1247945;
        let b = 3.1213820;
        let c = 2.5609104;
        let delta = 78.7590539;
        let mu = 0.6966326;
        let Qo = 312.1341346;
        let palp = 1.4074424;
        let bet = 0.0070975;
        let alp = 3.1083847;

        Self {
            A,
            B,
            rh,
            a,
            sig,
            lam,
            gam,
            b,
            c,
            delta,
            mu,
            Qo,
            palp,
            bet,
            alp,
            bg: a,
            eta: delta / Qo,
        }
    }
}
// f46cb3ad ends here

// [[file:../edip.note::3d96dcdb][3d96dcdb]]
pub use crate::edip::*;
// 3d96dcdb ends here

// [[file:../edip.note::4a762dc4][4a762dc4]]
#[cfg(feature = "adhoc")]
/// Docs for local mods
pub mod docs {
    macro_rules! export_doc {
        ($l:ident) => {
            pub mod $l {
                pub use crate::$l::*;
            }
        };
    }

    export_doc!(edip);
}
// 4a762dc4 ends here
