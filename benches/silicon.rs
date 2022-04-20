// [[file:../edip.note::f00f129e][f00f129e]]
use gchemol::prelude::*;
use gchemol::Molecule;
use gut::prelude::*;
use std::collections::{HashMap, HashSet};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use edip::*;

fn init() -> Result<(Vec<HashSet<usize>>, HashMap<(usize, usize), [f64; 3]>)> {
    let f = "./tests/files/silicon.xyz";
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

    Ok((neighbors, distances))
}

fn silicon(neighbors: &[HashSet<usize>], distances: &HashMap<(usize, usize), [f64; 3]>) {
    let n = neighbors.len();
    let mut f = vec![[0.0; 3]; n];
    let params = EdipParameters::silicon();
    let (_energy, _virial) = compute_forces(&mut f, &neighbors, &distances, &params);
}

fn criterion_benchmark(c: &mut Criterion) {
    let (neighbors, distances) = init().unwrap();
    c.bench_function("silicon", |b| b.iter(|| black_box(silicon(&neighbors, &distances))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
// f00f129e ends here
