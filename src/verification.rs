use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;

pub fn generate_phi() -> Vec<i64> {
    let mut rng = rand::thread_rng();
    let mut phi_k: Vec<i64> = Vec::new();
    for i in 0..L {
        phi_k.push(rng.gen_range(0..q));
    }
    phi_k
}
        

pub fn generate_omega() -> Vec<i64> {
    let mut rng = rand::thread_rng();
    let mut omega_k: Vec<i64> = Vec::new();
    for i in 0..256 {
        omega_k.push(rng.gen_range(0..q));
    }
    omega_k
}


pub fn sample_jl_projection() -> Array2<i64> {

    let between = Uniform::from(-1..=1);

    let mut rng = rand::thread_rng();

    let mut pi_i : Array2<i64> = Array2::zeros((256, N));

    for ((i, j), value) in pi_i.indexed_iter_mut() {
        *value = between.sample(&mut rng);
        //println!("matrix[{}][{}] = {}", i, j, value);
    }
    pi_i
}


pub fn valid_projection(projection: Array2<Polynomial<i64>>) -> bool {

    let val : f64 = 128.;
    let total_norm = compute_total_norm(projection);
    println!("TOTAL NORM OF JL PROJECTION: {}, {}",total_norm, val.sqrt()*(BETA_BOUND as f64));
    return total_norm <= val.sqrt()*(BETA_BOUND as f64);
}
