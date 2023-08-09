use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;


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


/*
pub fn valid_projection(projection: Array2<Polynomial<i64>>) : bool {

    compute_norm(&DVector::from_vec(projected_vec)) <= sqrt(128)*BETA_BOUND
}
*/
