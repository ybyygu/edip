// [[file:../edip.note::83a3f290][83a3f290]]
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

//! Literature
//!
//! http://www-math.mit.edu/~bazant/EDIP
//! M.Z. Bazant & E. Kaxiras: Modeling of Covalent Bonding in Solids by
//!                           Inversion of Cohesive Energy Curves;
//!                           Phys. Rev. Lett. 77, 4370 (1996)
//! M.Z. Bazant, E. Kaxiras and J.F. Justo: Environment-dependent interatomic
//!                                         potential for bulk silicon;
//!                                         Phys. Rev. B 56, 8542-8552 (1997)
// 83a3f290 ends here

// [[file:../edip.note::178e12ff][178e12ff]]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
use super::*;

use gut::prelude::*;
use std::collections::HashSet;

type Array3 = [f64; 3];
// 178e12ff ends here

// [[file:../edip.note::8a0e64b1][8a0e64b1]]
// tau(Z) (Ismail & Kaxiras, 1993)
const u1: f64 = -0.165799;
const u2: f64 = 32.557;
const u3: f64 = 0.286198;
const u4: f64 = 0.66;

/// Parameters for EDIP
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
// 8a0e64b1 ends here

// [[file:../edip.note::d9088664][d9088664]]
/* DEBUGGING FLAGS - systematically turn on various pieces of the potential. */
const V2_on: bool = true;
const V2Z_on: bool = true;
const V3_on: bool = true;
const V3g_on: bool = true;
const V3h_on: bool = true;
const V3Z_on: bool = true;

#[derive(Default, Clone)]
struct Store2 {
    // various V2 functions and derivatives
    t0: f64,
    t1: f64,
    t2: f64,
    t3: f64,
    // unit separation vector
    dx: f64,
    dy: f64,
    dz: f64,
    // bond length (only needed for virial)
    r: f64,
}

#[derive(Default, Clone)]
struct Store3 {
    // 3-body radial function and its derivative
    g: f64,
    dg: f64,
    // 1/r
    rinv: f64,
    // unit separation vector
    dx: f64,
    dy: f64,
    dz: f64,
    // bond length (only needed for virial)
    r: f64,
}

#[derive(Default, Clone)]
struct Storez {
    // derivative of neighbor function f'(r)
    df: f64,
    // array to accumulate coordination force prefactors
    sum: f64,
    // unit separation vector
    dx: f64,
    dy: f64,
    dz: f64,
    // bond length (only needed for virial)
    r: f64,
}
// d9088664 ends here

// [[file:../edip.note::9c60872a][9c60872a]]
fn compute_forces_edip(
    // Position of each atom.
    positions: &[Array3],
    // On output, contains the total force on each atom.
    forces: &mut [Array3],
    // A list of double counted neighbors of atom i (same indices as in the positions and forces ).
    neighbors: &[HashSet<usize>],
    // Parameter set for EDIP
    params: &EdipParams,
) -> (f64, f64) {
    // the total number of particles
    let n_own = positions.len();
    assert_eq!(forces.len(), n_own);
    assert_eq!(neighbors.len(), n_own);
    assert!(!neighbors.is_empty(), "invalid neighbors: {neighbors:?}");
    let max_nbrs = neighbors.iter().map(|x| x.len()).max().unwrap();

    /* atom ID numbers for s2[] */
    let mut num2 = vec![0; max_nbrs];
    /* atom ID numbers for s3[] */
    let mut num3 = vec![0; max_nbrs];
    /* atom ID numbers for sz[] */
    let mut numz = vec![0; max_nbrs];

    let mut s2 = vec![Store2::default(); max_nbrs];
    let mut s3 = vec![Store3::default(); max_nbrs];
    /* coordination number stuff, c<r<b */
    let mut sz = vec![Storez::default(); max_nbrs];

    /* INITIALIZE FORCES AND GLOBAL SUMS */
    for i in 0..n_own {
        forces[i] = [0.0; 3];
    }

    /* measurements in force routine */
    let mut e_potential = 0.0;
    // e_potential = e_v2 + e_v3
    let mut e_v2 = 0.0;
    let mut e_v3 = 0.0;
    let mut coord_total = 0.0;
    let mut virial = 0.0;

    /* COMBINE COEFFICIENTS */
    let asqr = params.a.powi(2);
    let Qort = params.Qo.sqrt();
    let muhalf = params.mu * 0.5;
    let u5 = u2 * u4;
    let bmc = params.b - params.c;
    let cmbinv = 1.0 / (params.c - params.b);

    /*--- LEVEL 1: OUTER LOOP OVER ATOMS ---*/
    for i in 0..n_own {
        /* RESET COORDINATION AND NEIGHBOR NUMBERS */
        let mut Z = 0.0;
        let mut n2 = 0; /* size of s2[] */
        let mut n3 = 0; /* size of s3[] */
        let mut nz = 0; /* size of sz[] */

        /*--- LEVEL 2: LOOP PREPASS OVER PAIRS ---*/
        for &j in &neighbors[i] {
            /* TEST IF WITHIN OUTER CUTOFF */
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            // dx = MIN_IMAGE_DISTANCE(dx,box.L_x,box.L_x_div2);
            // dy = MIN_IMAGE_DISTANCE(dy,box.L_y,box.L_y_div2);
            // dz = MIN_IMAGE_DISTANCE(dz,box.L_z,box.L_z_div2);
            let rsqr = dx * dx + dy * dy + dz * dz;
            if dx.abs() < params.a && dy.abs() < params.a && dz.abs() < params.a && rsqr < asqr {
                let r = rsqr.sqrt();
                /* PARTS OF TWO-BODY INTERACTION r<a */
                num2[n2] = j;
                let rinv = 1.0 / r;
                let dxrinv = dx * rinv;
                let dyrinv = dy * rinv;
                let dzrinv = dz * rinv;
                let rmainv = 1.0 / (r - params.a);
                s2[n2].t0 = params.A * (params.sig * rmainv).exp();
                s2[n2].t1 = (params.B * rinv).powf(params.rh);
                s2[n2].t2 = params.rh * rinv;
                s2[n2].t3 = params.sig * rmainv * rmainv;
                s2[n2].dx = dxrinv;
                s2[n2].dy = dyrinv;
                s2[n2].dz = dzrinv;
                s2[n2].r = r;
                n2 += 1;

                /* RADIAL PARTS OF THREE-BODY INTERACTION r<b */
                if r < params.bg {
                    num3[n3] = j;
                    let rmbinv = 1.0 / (r - params.bg);
                    let temp1 = params.gam * rmbinv;
                    let temp0 = temp1.exp();
                    if V3g_on {
                        s3[n3].g = temp0;
                        s3[n3].dg = -rmbinv * temp1 * temp0;
                    } else {
                        s3[n3].g = 1.0;
                        s3[n3].dg = 0.0;
                    }
                    s3[n3].dx = dxrinv;
                    s3[n3].dy = dyrinv;
                    s3[n3].dz = dzrinv;
                    s3[n3].rinv = rinv;
                    s3[n3].r = r;
                    n3 += 1;

                    /* COORDINATION AND NEIGHBOR FUNCTION c<r<b */
                    if r < params.b {
                        if r < params.c {
                            Z += 1.0;
                        } else {
                            let xinv = bmc / (r - params.c);
                            let xinv3 = xinv.powi(3);
                            let den = 1.0 / (1.0 - xinv3);
                            let temp1 = params.alp * den;
                            let fZ = temp1.exp();
                            Z += fZ;
                            numz[nz] = j;
                            sz[nz].df = fZ * temp1 * den * 3.0 * xinv3 * xinv * cmbinv; /* df/dr */
                            sz[nz].dx = dxrinv;
                            sz[nz].dy = dyrinv;
                            sz[nz].dz = dzrinv;
                            sz[nz].r = r;
                            nz += 1;
                        }
                    }
                }
            }
        }
        coord_total += Z;

        /* ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES */
        for nl in 0..nz {
            sz[nl].sum = 0.0;
        }

        /* ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION */
        let mut pz = 0.0;
        let mut dp = 0.0;
        if V2_on {
            if V2Z_on {
                let temp0 = params.bet * Z;
                pz = params.palp * (-temp0 * Z).exp(); /* bond order */
                dp = -2.0 * temp0 * pz; /* derivative of bond order */
            } else {
                pz = params.palp * (-params.bet * 16.0).exp();
            }
            /*--- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---*/
            for nj in 0..n2 {
                let temp0 = s2[nj].t1 - pz;
                /* two-body energy V2(rij,Z) */
                e_v2 += temp0 * s2[nj].t0;

                /* two-body forces */
                let dV2j = -(s2[nj].t0) * ((s2[nj].t1) * (s2[nj].t2) + temp0 * (s2[nj].t3)); /* dV2/dr */
                let dV2ijx = dV2j * s2[nj].dx;
                let dV2ijy = dV2j * s2[nj].dy;
                let dV2ijz = dV2j * s2[nj].dz;
                forces[i][0] += dV2ijx;
                forces[i][1] += dV2ijy;
                forces[i][2] += dV2ijz;
                let j = num2[nj];
                forces[j][0] -= dV2ijx;
                forces[j][1] -= dV2ijy;
                forces[j][2] -= dV2ijz;

                /* dV2/dr contribution to virial */
                virial -= s2[nj].r * (dV2ijx * s2[nj].dx + dV2ijy * s2[nj].dy + dV2ijz * s2[nj].dz);

                /*--- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---*/
                let dV2dZ = -dp * s2[nj].t0;
                for nl in 0..nz {
                    sz[nl].sum += dV2dZ;
                }
            }
        }

        /* COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION */
        let mut tau = 0.0;
        let mut dtau = 0.0;
        let mut winv = 0.0;
        let mut dwinv = 0.0;
        if V3_on {
            if V3Z_on {
                winv = Qort * (-muhalf * Z).exp(); /* inverse width of angular function */
                dwinv = -muhalf * winv; /* its derivative */
                let temp0 = (-u4 * Z).exp();
                tau = u1 + u2 * temp0 * (u3 - temp0); /* -cosine of angular minimum */
                dtau = u5 * temp0 * (2.0 * temp0 - u3); /* its derivative */
            } else {
                winv = Qort * (-muhalf * 4.0).exp();
                dwinv = 0.0;
                tau = 1.0 / 3.0;
                dtau = 0.0;
            }

            /*--- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---*/
            for nj in 0..(n3 - 1) {
                let j = num3[nj];
                /*--- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---*/
                for nk in (nj + 1)..n3 {
                    let k = num3[nk];
                    /* angular function h(l,Z) */
                    let lcos = s3[nj].dx * s3[nk].dx + s3[nj].dy * s3[nk].dy + s3[nj].dz * s3[nk].dz;
                    let x = (lcos + tau) * winv;
                    let temp0 = (-x * x).exp();

                    let H;
                    let dhdl;
                    let mut dHdx = 0.0;
                    if V3h_on {
                        H = params.lam * (1.0 - temp0 + params.eta * x * x);
                        dHdx = 2.0 * params.lam * x * (temp0 + params.eta);
                        dhdl = dHdx * winv;
                    } else {
                        H = 1.0;
                        dhdl = 0.0;
                    }

                    /* three-body energy */
                    let temp1 = s3[nj].g * s3[nk].g;
                    e_v3 += temp1 * H;

                    /* (-) radial force on atom j */
                    let dV3rij = s3[nj].dg * s3[nk].g * H;
                    let dV3rijx = dV3rij * s3[nj].dx;
                    let dV3rijy = dV3rij * s3[nj].dy;
                    let dV3rijz = dV3rij * s3[nj].dz;
                    let mut fjx = dV3rijx;
                    let mut fjy = dV3rijy;
                    let mut fjz = dV3rijz;

                    /* (-) radial force on atom k */
                    let dV3rik = s3[nj].g * s3[nk].dg * H;
                    let dV3rikx = dV3rik * s3[nk].dx;
                    let dV3riky = dV3rik * s3[nk].dy;
                    let dV3rikz = dV3rik * s3[nk].dz;
                    let mut fkx = dV3rikx;
                    let mut fky = dV3riky;
                    let mut fkz = dV3rikz;

                    /* (-) angular force on j */
                    let dV3l = temp1 * dhdl;
                    let dV3ljx = dV3l * (s3[nk].dx - lcos * s3[nj].dx) * s3[nj].rinv;
                    let dV3ljy = dV3l * (s3[nk].dy - lcos * s3[nj].dy) * s3[nj].rinv;
                    let dV3ljz = dV3l * (s3[nk].dz - lcos * s3[nj].dz) * s3[nj].rinv;
                    fjx += dV3ljx;
                    fjy += dV3ljy;
                    fjz += dV3ljz;

                    /* (-) angular force on k */
                    let dV3lkx = dV3l * (s3[nj].dx - lcos * s3[nk].dx) * s3[nk].rinv;
                    let dV3lky = dV3l * (s3[nj].dy - lcos * s3[nk].dy) * s3[nk].rinv;
                    let dV3lkz = dV3l * (s3[nj].dz - lcos * s3[nk].dz) * s3[nk].rinv;
                    fkx += dV3lkx;
                    fky += dV3lky;
                    fkz += dV3lkz;

                    /* apply radial + angular forces to i, j, k */
                    forces[j][0] -= fjx;
                    forces[j][1] -= fjy;
                    forces[j][2] -= fjz;
                    forces[k][0] -= fkx;
                    forces[k][1] -= fky;
                    forces[k][2] -= fkz;
                    forces[i][0] += fjx + fkx;
                    forces[i][1] += fjy + fky;
                    forces[i][2] += fjz + fkz;

                    /* dV3/dR contributions to virial */
                    virial -= s3[nj].r * (fjx * s3[nj].dx + fjy * s3[nj].dy + fjz * s3[nj].dz);
                    virial -= s3[nk].r * (fkx * s3[nk].dx + fky * s3[nk].dy + fkz * s3[nk].dz);

                    /* prefactor for 4-body forces from coordination */
                    if V3Z_on {
                        let dxdZ = dwinv * (lcos + tau) + winv * dtau;
                        let dV3dZ = temp1 * dHdx * dxdZ;
                        /*--- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---*/
                        for nl in 0..nz {
                            sz[nl].sum += dV3dZ;
                        }
                    }
                }
            }
        } // end if V3_on

        /*--- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---*/
        for nl in 0..nz {
            let dedrl = sz[nl].sum * sz[nl].df;
            let dedrlx = dedrl * sz[nl].dx;
            let dedrly = dedrl * sz[nl].dy;
            let dedrlz = dedrl * sz[nl].dz;
            forces[i][0] += dedrlx;
            forces[i][1] += dedrly;
            forces[i][2] += dedrlz;
            let l = numz[nl];
            forces[l][0] -= dedrlx;
            forces[l][1] -= dedrly;
            forces[l][2] -= dedrlz;

            /* dE/dZ*dZ/dr contribution to virial */
            virial -= sz[nl].r * (dedrlx * sz[nl].dx + dedrly * sz[nl].dy + dedrlz * sz[nl].dz);
        }
    }
    let e_potential = e_v2 + e_v3;
    virial /= 3.0;

    (e_potential, virial)
}
// 9c60872a ends here

// [[file:../edip.note::bccd1cf0][bccd1cf0]]
#[test]
fn test_edip() -> Result<()> {
    use gchemol::prelude::*;
    use gchemol::Molecule;
    use vecfx::*;

    let f = "./tests/files/si5.xyz";
    let mol = Molecule::from_file(f)?;

    let n = mol.natoms();
    let pos = mol.positions().collect_vec();
    let mut neighbors: Vec<HashSet<usize>> = vec![];
    for i in 0..n {
        let mut xx: HashSet<_> = (0..n).collect();
        xx.remove(&i);
        neighbors.push(xx.into_iter().collect());
    }
    let mut f = vec![[0.0; 3]; n];
    let params = EdipParams::silicon();
    let (energy, virial) = compute_forces_edip(&pos, &mut f, &neighbors, &params);
    approx::assert_relative_eq!(energy, -14.566606, epsilon = 1e-5);
    approx::assert_relative_eq!(virial, -3.643552, epsilon = 1e-5);
    #[rustfmt::skip]
    let f_expected = [ -0.19701000,  -0.62522600,   0.02948000,
                       -0.23330600,   0.43698600,   0.42511400,
                        0.43297600,   0.18333500,  -0.44864300,
                       -1.75669300,   0.50149400,  -1.48879300,
                        1.75403300,  -0.49658900,   1.48284200];
    approx::assert_relative_eq!(
        f.as_flat().as_vector_slice(),
        f_expected.as_vector_slice(),
        epsilon = 1E-5,
    );

    Ok(())
}
// bccd1cf0 ends here
