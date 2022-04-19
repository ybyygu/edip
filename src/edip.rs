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

use std::collections::{HashMap, HashSet};
// 178e12ff ends here

// [[file:../edip.note::d9088664][d9088664]]
#[derive(Default, Clone, Copy)]
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

#[derive(Default, Clone, Copy)]
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

#[derive(Default, Clone, Copy)]
struct Storez {
    // derivative of neighbor function f'(r)
    df: f64,
    // unit separation vector
    dx: f64,
    dy: f64,
    dz: f64,
    // bond length (only needed for virial)
    r: f64,
}
// d9088664 ends here

// [[file:../edip.note::44c1697e][44c1697e]]
impl EdipParams {
    fn compute_prepass(
        self,
        s2: &mut HashMap<usize, Store2>,
        s3: &mut HashMap<usize, Store3>,
        sz: &mut HashMap<usize, Storez>,
        neighbors: &[HashSet<usize>],
        distances: &HashMap<(usize, usize), Array3>,
        i: usize,
    ) -> f64 {
        // Reset Coordination And Neighbor Numbers
        let mut Z = 0.0;

        let cmb_inv = 1.0 / (self.c - self.b);
        let asqr = self.a.powi(2);
        s2.clear();
        s3.clear();
        sz.clear();

        // Get all pairs within cutoff in parallel
        let pairs_within_cutoff: Vec<_> = neighbors[i]
            .par_iter()
            .filter_map(|&j| {
                assert_ne!(i, j);
                let [dx, dy, dz] = distances[&(i, j)];
                if dx.abs() < self.a && dy.abs() < self.a && dz.abs() < self.a {
                    let rsqr = dx * dx + dy * dy + dz * dz;
                    if rsqr < asqr {
                        let r = rsqr.sqrt();
                        return (i, j, r, dx, dy, dz).into();
                    }
                }
                None
            })
            .collect();

        // parts of two-body interaction r<a
        for &(i, j, r, dx, dy, dz) in pairs_within_cutoff.iter() {
            let rinv = 1.0 / r;
            let dxrinv = dx * rinv;
            let dyrinv = dy * rinv;
            let dzrinv = dz * rinv;
            let rmainv = 1.0 / (r - self.a);
            let temp = Store2 {
                t0: self.A * (self.sig * rmainv).exp(),
                t1: (self.B * rinv).powf(self.rh),
                t2: self.rh * rinv,
                t3: self.sig * rmainv * rmainv,
                dx: dxrinv,
                dy: dyrinv,
                dz: dzrinv,
                r,
            };
            s2.insert(j, temp);
        }

        // radial parts of three-body interaction r<b
        for (&j, &Store2 { dx, dy, dz, r, .. }) in s2.iter() {
            if r < self.bg {
                let rinv = 1.0 / r;
                let rmbinv = 1.0 / (r - self.bg);
                let temp1 = self.gam * rmbinv;
                let temp0 = temp1.exp();
                let temp = Store3 {
                    dx,
                    dy,
                    dz,
                    dg: -rmbinv * temp1 * temp0,
                    g: temp0,
                    rinv,
                    r,
                };
                s3.insert(j, temp);
            }
        }

        // coordination and neighbor function c<r<b
        for (&j, &Store3 { dx, dy, dz, r, .. }) in s3.iter() {
            if r < self.b {
                if r < self.c {
                    Z += 1.0;
                } else {
                    let xinv = (self.b - self.c) / (r - self.c);
                    let xinv3 = xinv.powi(3);
                    let den = 1.0 / (1.0 - xinv3);
                    let temp1 = self.alp * den;
                    let fZ = temp1.exp();
                    Z += fZ;
                    let temp = Storez {
                        dx,
                        dy,
                        dz,
                        df: fZ * temp1 * den * 3.0 * xinv3 * xinv * cmb_inv, /* df/dr */
                        r,
                    };
                    sz.insert(j, temp);
                }
            }
        }

        Z
    }
}
// 44c1697e ends here

// [[file:../edip.note::0c7bb6bc][0c7bb6bc]]
/// Environment-dependence of pair interaction
struct PairInteraction {
    /// bond order
    pz: f64,
    /// derivative of bond order
    dp: f64,
}

impl PairInteraction {
    fn compute(&self, s2: &HashMap<usize, Store2>, forces: &mut [Array3], i: usize) -> (f64, f64, f64) {
        let pz = self.pz;
        let dp = self.dp;

        let mut virial = 0.0;
        for (&j, s2nj) in s2 {
            // two-body forces
            // dV2/dr
            let dV2j = -s2nj.t0 * (s2nj.t1 * s2nj.t2 + (s2nj.t1 - pz) * s2nj.t3);
            let dV2ijx = dV2j * s2nj.dx;
            let dV2ijy = dV2j * s2nj.dy;
            let dV2ijz = dV2j * s2nj.dz;
            forces[i][0] += dV2ijx;
            forces[i][1] += dV2ijy;
            forces[i][2] += dV2ijz;
            forces[j][0] -= dV2ijx;
            forces[j][1] -= dV2ijy;
            forces[j][2] -= dV2ijz;

            // dV2/dr contribution to virial
            virial -= s2nj.r * (dV2ijx * s2nj.dx + dV2ijy * s2nj.dy + dV2ijz * s2nj.dz);
        }

        // for pair coordination force prefactors
        let sum = s2.values().map(|x| -dp * x.t0).sum();
        // two-body energy V2(rij,Z)
        let energy = s2.values().map(|x| (x.t1 - pz) * x.t0).sum();

        (energy, virial, sum)
    }
}
// 0c7bb6bc ends here

// [[file:../edip.note::1c360840][1c360840]]
struct ThreeBodyInteraction {
    tau: f64,
    w_inv: f64,
    dw_inv: f64,
    dtau: f64,
    temp0: f64,
}

impl ThreeBodyInteraction {
    fn compute(
        &self,
        s3: &HashMap<usize, Store3>,
        forces: &mut [Array3],
        i: usize,
        params: &EdipParams,
    ) -> (f64, f64, f64) {
        let tau = self.tau;
        let w_inv = self.w_inv;
        let dw_inv = self.dw_inv;
        let dtau = self.dtau;
        let temp0 = self.temp0;

        // --- Level 2: First Loop For Three-Body Interactions ---
        let mut virial = 0.0;
        let mut e_v3 = 0.0;
        // accumulate coordination force prefactors
        let mut sum = 0.0;
        for p in s3.iter().combinations(2) {
            let (&j, s3nj) = p[0];
            let (&k, s3nk) = p[1];
            /* angular function h(l,Z) */
            let lcos = s3nj.dx * s3nk.dx + s3nj.dy * s3nk.dy + s3nj.dz * s3nk.dz;
            let x = (lcos + tau) * w_inv;
            let temp0 = (-x * x).exp();

            let H = params.lam * (1.0 - temp0 + params.eta * x * x);
            let dHdx = 2.0 * params.lam * x * (temp0 + params.eta);
            let dhdl = dHdx * w_inv;

            // three-body energy
            let temp1 = s3nj.g * s3nk.g;
            e_v3 += temp1 * H;

            // (-) radial force on atom j
            let dV3rij = s3nj.dg * s3nk.g * H;
            let dV3rijx = dV3rij * s3nj.dx;
            let dV3rijy = dV3rij * s3nj.dy;
            let dV3rijz = dV3rij * s3nj.dz;
            let mut fjx = dV3rijx;
            let mut fjy = dV3rijy;
            let mut fjz = dV3rijz;

            // (-) radial force on atom k
            let dV3rik = s3nj.g * s3nk.dg * H;
            let dV3rikx = dV3rik * s3nk.dx;
            let dV3riky = dV3rik * s3nk.dy;
            let dV3rikz = dV3rik * s3nk.dz;
            let mut fkx = dV3rikx;
            let mut fky = dV3riky;
            let mut fkz = dV3rikz;

            // (-) angular force on j
            let dV3l = temp1 * dhdl;
            let dV3ljx = dV3l * (s3nk.dx - lcos * s3nj.dx) * s3nj.rinv;
            let dV3ljy = dV3l * (s3nk.dy - lcos * s3nj.dy) * s3nj.rinv;
            let dV3ljz = dV3l * (s3nk.dz - lcos * s3nj.dz) * s3nj.rinv;
            fjx += dV3ljx;
            fjy += dV3ljy;
            fjz += dV3ljz;

            // (-) angular force on k
            let dV3lkx = dV3l * (s3nj.dx - lcos * s3nk.dx) * s3nk.rinv;
            let dV3lky = dV3l * (s3nj.dy - lcos * s3nk.dy) * s3nk.rinv;
            let dV3lkz = dV3l * (s3nj.dz - lcos * s3nk.dz) * s3nk.rinv;
            fkx += dV3lkx;
            fky += dV3lky;
            fkz += dV3lkz;

            // apply radial + angular forces to i, j, k
            forces[j][0] -= fjx;
            forces[j][1] -= fjy;
            forces[j][2] -= fjz;
            forces[k][0] -= fkx;
            forces[k][1] -= fky;
            forces[k][2] -= fkz;
            forces[i][0] += fjx + fkx;
            forces[i][1] += fjy + fky;
            forces[i][2] += fjz + fkz;

            // dV3/dR contributions to virial
            virial -= s3nj.r * (fjx * s3nj.dx + fjy * s3nj.dy + fjz * s3nj.dz);
            virial -= s3nk.r * (fkx * s3nk.dx + fky * s3nk.dy + fkz * s3nk.dz);

            // Prefactor for 4-body forces from coordination
            let dxdZ = dw_inv * (lcos + tau) + w_inv * dtau;
            let dV3dZ = temp1 * dHdx * dxdZ;
            // Three-Body Coordination Forces
            sum += dV3dZ;
        }

        (e_v3, virial, sum)
    }
}
// 1c360840 ends here

// [[file:../edip.note::9c60872a][9c60872a]]
pub fn compute_forces_edip(
    // The total forces to be computed
    forces: &mut [Array3],
    // A list of double counted neighbors of atom i (same indices as in the forces ).
    neighbors: &[HashSet<usize>],
    // Contains relative coordinates of atom j relative to its neighboring atom
    // i. Will panic if no data for any pair of neighboring atoms i and j
    distances: &HashMap<(usize, usize), Array3>,
    // Parameter set for EDIP
    params: &EdipParams,
) -> (f64, f64) {
    // the total number of particles
    let n_own = forces.len();
    assert_eq!(neighbors.len(), n_own);
    assert!(!neighbors.is_empty(), "invalid neighbors: {neighbors:?}");
    let max_nbrs = neighbors.iter().map(|x| x.len()).max().unwrap();

    let mut s2 = HashMap::new();
    let mut s3 = HashMap::new();
    // Coordination number stuff, c<r<b
    let mut sz = HashMap::new();

    // Initialize forces
    for i in 0..n_own {
        forces[i] = [0.0; 3];
    }

    // e_potential = energy_v2 + energy_v3
    let mut energy_v2 = 0.0;
    let mut energy_v3 = 0.0;
    let mut virial = 0.0;

    // outer loop over atoms
    for i in 0..n_own {
        // reset coordination and neighbor numbers
        let Z = params.compute_prepass(&mut s2, &mut s3, &mut sz, neighbors, distances, i);

        // environment-dependence of pair interaction
        let temp0 = params.bet * Z;
        // bond order
        let pz = params.palp * (-temp0 * Z).exp();
        // derivative of bond order
        let dp = -2.0 * temp0 * pz;
        let pair_int = PairInteraction { pz, dp };
        let (e, v, sum_s2) = pair_int.compute(&s2, forces, i);
        energy_v2 += e;
        virial += v;

        // coordination-dependence of three-body interaction
        let temp0 = (-u4 * Z).exp();
        let w_inv = params.Qo.sqrt() * (-0.5 * params.mu * Z).exp(); /* inverse width of angular function */
        let threebody = ThreeBodyInteraction {
            temp0,
            w_inv,
            tau: u1 + u2 * temp0 * (u3 - temp0), /* -cosine of angular minimum */
            dtau: u2 * u4 * temp0 * (2.0 * temp0 - u3), /* its derivative */
            dw_inv: -0.5 * params.mu * w_inv,    /* its derivative */
        };
        let (e, v, sum_s3) = threebody.compute(&s3, forces, i, &params);
        energy_v3 += e;
        virial += v;

        // loop to apply coordination forces
        let psum = sum_s2 + sum_s3;
        for (&l, sznl) in &sz {
            let dedrl = psum * sznl.df;
            let dedrlx = dedrl * sznl.dx;
            let dedrly = dedrl * sznl.dy;
            let dedrlz = dedrl * sznl.dz;
            forces[i][0] += dedrlx;
            forces[i][1] += dedrly;
            forces[i][2] += dedrlz;
            forces[l][0] -= dedrlx;
            forces[l][1] -= dedrly;
            forces[l][2] -= dedrlz;

            // dE/dZ*dZ/dr contribution to virial
            virial -= sznl.r * (dedrlx * sznl.dx + dedrly * sznl.dy + dedrlz * sznl.dz);
        }
    }
    let e_potential = energy_v2 + energy_v3;
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
    let mut distances = HashMap::new();
    for i in 0..n {
        let mut xx: HashSet<_> = (0..n).collect();
        xx.remove(&i);
        neighbors.push(xx.into_iter().collect());
        for j in 0..n {
            if i != j {
                let pi = positions[i];
                let pj = positions[j];
                let dx = pj[0] - pi[0];
                let dy = pj[1] - pi[1];
                let dz = pj[2] - pi[2];
                distances.insert((i, j), [dx, dy, dz]);
            }
        }
    }

    let mut f = vec![[0.0; 3]; n];
    let params = EdipParams::silicon();
    let (energy, virial) = compute_forces_edip(&mut f, &neighbors, &distances, &params);
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
