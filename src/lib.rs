// [[file:../edip.note::c5791ec5][c5791ec5]]
//! # Literature
//!
//! - <http://www-math.mit.edu/~bazant/EDIP>
//!
//! - M.Z. Bazant & E. Kaxiras: Modeling of Covalent Bonding in Solids by
//! Inversion of Cohesive Energy Curves; Phys. Rev. Lett. 77, 4370 (1996)
//!
//! - M.Z. Bazant, E. Kaxiras and J.F. Justo: Environment-dependent interatomic
//! potential for bulk silicon; Phys. Rev. B 56, 8542-8552 (1997)
// c5791ec5 ends here

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
/* EDIP Si PARAMETERS Justo et al., Phys. Rev. B 58, 2539 (1998).

     5.6714030     2.0002804     1.2085196     3.1213820     0.5774108
     1.4533108     1.1247945     3.1213820     2.5609104    78.7590539
     0.6966326   312.1341346     1.4074424     0.0070975     3.1083847

connection between these parameters and Justo et al., Phys. Rev. B 58, 2539 (1998):

A((B/r)**rh-palp*exp(-bet*Z*Z)) = A'((B'/r)**rh-exp(-bet*Z*Z))

so in the paper (')
A' = A*palp
B' = B * palp**(-1/rh)
eta = detla/Qo
*/

// tau(Z) (Ismail & Kaxiras, 1993)
const u1: f64 = -0.165799;
const u2: f64 = 32.557;
const u3: f64 = 0.286198;
const u4: f64 = 0.66;

/// Parameters for EDIP
#[derive(Debug, Copy, Clone)]
pub struct EdipParameters {
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
    eta: f64,
}

impl EdipParameters {
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
