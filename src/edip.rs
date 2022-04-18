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
use super::*;

use gut::prelude::*;
use std::collections::HashSet;

const MAX_PART: usize = 4097; /* maximum number of particles (no dynamic mem.alloc.) */
const MAX_NBRS: usize = 30;

#[derive(Debug, Clone, Default)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}
// 178e12ff ends here

// [[file:../edip.note::8a0e64b1][8a0e64b1]]
// tau(Z) (Ismail & Kaxiras, 1993)
const u1: f64 = -0.165799;
const u2: f64 = 32.557;
const u3: f64 = 0.286198;
const u4: f64 = 0.66;

struct EdipParams {
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
    fn silicon() -> Self {
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
const Zfast: bool = true;

#[derive(Debug, Default, Clone)]
struct store2 {
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

#[derive(Debug, Default, Clone)]
struct store3 {
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

#[derive(Debug, Default, Clone)]
struct storez {
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

fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

fn fabs(x: f64) -> f64 {
    x.abs()
}

fn exp(x: f64) -> f64 {
    x.exp()
}

fn pow(x: f64, r: f64) -> f64 {
    x.powf(r)
}
// d9088664 ends here

// [[file:../edip.note::9c60872a][9c60872a]]
fn compute_forces_edip(
    // Position of each atom.
    pos: &[Vector3],
    // On output, contains the total force on each atom.
    f: &mut [Vector3],
    // A list of neighboring atoms of atom i (same indices as in the pos and f
    // ).
    neighbors: &[HashSet<usize>],
) -> (f64, f64) {
    /* measurements in force routine */
    let mut virial = 0.0;
    let mut e_potential = 0.0;
    let mut v2 = 0.0;
    let mut v3 = 0.0;
    let mut v2sum = 0.0;
    let mut coord_total = 0.0;

    /* atom ID numbers for s2[] */
    let mut num2 = vec![0; MAX_NBRS];
    /* atom ID numbers for s3[] */
    let mut num3 = vec![0; MAX_NBRS];
    /* atom ID numbers for sz[] */
    let mut numz = vec![0; MAX_NBRS];

    let mut s2 = vec![store2::default(); MAX_NBRS];
    let mut s3 = vec![store3::default(); MAX_NBRS];
    /* coordination number stuff, c<r<b */
    let mut sz = vec![storez::default(); MAX_NBRS];
    let EdipParams {
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
        bg,
        eta,
    } = EdipParams::silicon();

    /* INITIALIZE FORCES AND GLOBAL SUMS */
    assert_eq!(pos.len(), f.len());
    // the total number of particles
    let N_own = pos.len();
    for i in 0..N_own {
        f[i].x = 0.0;
        f[i].y = 0.0;
        f[i].z = 0.0;
    }

    coord_total = 0.0;
    virial = 0.0;
    v2 = 0.0; /* e_potential = v2 + v3 */
    v3 = 0.0;

    /* COMBINE COEFFICIENTS */
    let asqr = a * a;
    let Qort = sqrt(Qo);
    let muhalf = mu * 0.5;
    let u5 = u2 * u4;
    let bmc = b - c;
    let cmbinv = 1.0 / (c - b);

    /*--- LEVEL 1: OUTER LOOP OVER ATOMS ---*/
    for i in 0..N_own {
        /* RESET COORDINATION AND NEIGHBOR NUMBERS */
        let mut Z = 0.0;
        let mut n2 = 0; /* size of s2[] */
        let mut n3 = 0; /* size of s3[] */
        let mut nz = 0; /* size of sz[] */

        /*--- LEVEL 2: LOOP PREPASS OVER PAIRS ---*/
        for &j in &neighbors[i] {
            /* TEST IF WITHIN OUTER CUTOFF */
            let mut dx = pos[j].x - pos[i].x;
            // dx = MIN_IMAGE_DISTANCE(dx,box.L_x,box.L_x_div2);
            if fabs(dx) < a {
                let mut dy = pos[j].y - pos[i].y;
                // dy = MIN_IMAGE_DISTANCE(dy,box.L_y,box.L_y_div2);
                if fabs(dy) < a {
                    let mut dz = pos[j].z - pos[i].z;
                    // dz = MIN_IMAGE_DISTANCE(dz,box.L_z,box.L_z_div2);
                    if fabs(dz) < a {
                        let rsqr = dx * dx + dy * dy + dz * dz;
                        if rsqr < asqr {
                            let r = sqrt(rsqr);
                            /* PARTS OF TWO-BODY INTERACTION r<a */
                            num2[n2] = j;
                            let rinv = 1.0 / r;
                            dx *= rinv;
                            dy *= rinv;
                            dz *= rinv;
                            let rmainv = 1.0 / (r - a);
                            s2[n2].t0 = A * exp(sig * rmainv);
                            s2[n2].t1 = pow(B * rinv, rh);
                            s2[n2].t2 = rh * rinv;
                            s2[n2].t3 = sig * rmainv * rmainv;
                            s2[n2].dx = dx;
                            s2[n2].dy = dy;
                            s2[n2].dz = dz;
                            s2[n2].r = r;
                            n2 += 1;

                            /* RADIAL PARTS OF THREE-BODY INTERACTION r<b */
                            if r < bg {
                                num3[n3] = j;
                                let rmbinv = 1.0 / (r - bg);
                                let temp1 = gam * rmbinv;
                                let temp0 = exp(temp1);
                                if V3g_on {
                                    s3[n3].g = temp0;
                                    s3[n3].dg = -rmbinv * temp1 * temp0;
                                } else {
                                    s3[n3].g = 1.0;
                                    s3[n3].dg = 0.0;
                                }
                                s3[n3].dx = dx;
                                s3[n3].dy = dy;
                                s3[n3].dz = dz;
                                s3[n3].rinv = rinv;
                                s3[n3].r = r;
                                n3 += 1;

                                /* COORDINATION AND NEIGHBOR FUNCTION c<r<b */
                                if r < b {
                                    if r < c {
                                        Z += 1.0;
                                    } else {
                                        let xinv = bmc / (r - c);
                                        let xinv3 = xinv * xinv * xinv;
                                        let den = 1.0 / (1.0 - xinv3);
                                        let temp1 = alp * den;
                                        let fZ = exp(temp1);
                                        Z += fZ;
                                        numz[nz] = j;
                                        sz[nz].df = fZ * temp1 * den * 3.0 * xinv3 * xinv * cmbinv; /* df/dr */
                                        sz[nz].dx = dx;
                                        sz[nz].dy = dy;
                                        sz[nz].dz = dz;
                                        sz[nz].r = r;
                                        nz += 1;
                                    }
                                }
                            }
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
            if V2_on {
                if V2Z_on {
                    let temp0 = bet * Z;
                    pz = palp * exp(-temp0 * Z); /* bond order */
                    dp = -2.0 * temp0 * pz; /* derivative of bond order */
                } else {
                    pz = palp * exp(-bet * 16.0);
                }
            }

            /*--- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---*/
            for nj in 0..n2 {
                let temp0 = s2[nj].t1 - pz;
                /* two-body energy V2(rij,Z) */
                v2 += temp0 * s2[nj].t0;

                /* two-body forces */
                let dV2j = -(s2[nj].t0) * ((s2[nj].t1) * (s2[nj].t2) + temp0 * (s2[nj].t3)); /* dV2/dr */
                let dV2ijx = dV2j * s2[nj].dx;
                let dV2ijy = dV2j * s2[nj].dy;
                let dV2ijz = dV2j * s2[nj].dz;
                f[i].x += dV2ijx;
                f[i].y += dV2ijy;
                f[i].z += dV2ijz;
                let j = num2[nj];
                f[j].x -= dV2ijx;
                f[j].y -= dV2ijy;
                f[j].z -= dV2ijz;

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
                winv = Qort * exp(-muhalf * Z); /* inverse width of angular function */
                dwinv = -muhalf * winv; /* its derivative */
                let temp0 = exp(-u4 * Z);
                tau = u1 + u2 * temp0 * (u3 - temp0); /* -cosine of angular minimum */
                dtau = u5 * temp0 * (2.0 * temp0 - u3); /* its derivative */
            } else {
                winv = Qort * exp(-muhalf * 4.0);
                dwinv = 0.0;
                tau = 1.0 / 3.0;
                dtau = 0.0;
            }

            /*--- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---*/
            for nj in 0..n3 - 1 {
                let j = num3[nj];
                /*--- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---*/
                for nk in (nj + 1)..n3 {
                    let k = num3[nk];
                    /* angular function h(l,Z) */
                    let lcos = s3[nj].dx * s3[nk].dx + s3[nj].dy * s3[nk].dy + s3[nj].dz * s3[nk].dz;
                    let x = (lcos + tau) * winv;
                    let temp0 = exp(-x * x);

                    let mut H;
                    let mut dhdl;
                    let mut dHdx = 0.0;
                    if V3h_on {
                        H = lam * (1.0 - temp0 + eta * x * x);
                        dHdx = 2.0 * lam * x * (temp0 + eta);
                        dhdl = dHdx * winv;
                    } else {
                        H = 1.0;
                        dhdl = 0.0;
                    }

                    /* three-body energy */
                    let temp1 = s3[nj].g * s3[nk].g;
                    v3 += temp1 * H;

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
                    f[j].x -= fjx;
                    f[j].y -= fjy;
                    f[j].z -= fjz;
                    f[k].x -= fkx;
                    f[k].y -= fky;
                    f[k].z -= fkz;
                    f[i].x += fjx + fkx;
                    f[i].y += fjy + fky;
                    f[i].z += fjz + fkz;

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
            f[i].x += dedrlx;
            f[i].y += dedrly;
            f[i].z += dedrlz;
            let l = numz[nl];
            f[l].x -= dedrlx;
            f[l].y -= dedrly;
            f[l].z -= dedrlz;

            /* dE/dZ*dZ/dr contribution to virial */
            virial -= sz[nl].r * (dedrlx * sz[nl].dx + dedrly * sz[nl].dy + dedrlz * sz[nl].dz);
        }
    }
    let e_potential = v2 + v3;
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
    let positions = mol.positions().collect_vec();
    let mut neighbors: Vec<HashSet<usize>> = vec![];
    for i in 0..n {
        let mut xx: HashSet<_> = (0..n).collect();
        xx.remove(&i);
        neighbors.push(xx.into_iter().collect());
    }
    let pos = mol.positions().map(|[x, y, z]| Vector3 { x, y, z }).collect_vec();
    let mut f = vec![Vector3 { x: 0.0, y: 0.0, z: 0.0 }; n];
    let (energy, virial) = compute_forces_edip(&pos, &mut f, &neighbors);
    approx::assert_relative_eq!(energy, -14.566606, epsilon = 1e-5);
    approx::assert_relative_eq!(virial, -3.643552, epsilon = 1e-5);

    let f = f.iter().map(|a| [a.x, a.y, a.z]).collect_vec();
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
