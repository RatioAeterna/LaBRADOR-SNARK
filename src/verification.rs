use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;

pub struct Verifier {

    projection : Array2<Polynomial<i64>>,
    phi_k : Vec<i64>,
    omega_k : Vec<i64>,
    b_prime : i64,

    // JL projection matrix
    Pi : Vec<Array2<i64>>,
}

impl Verifier {

    pub fn new() {}

    pub fn generate_phi() -> Vec<i64> {
        let mut rng = rand::thread_rng();
        let mut phi_k: Vec<i64> = Vec::new();
        for i in 0..L {
            phi_k.push(rng.gen_range(0..q));
        }
        phi_k
    }
            

    pub fn generate_omega(&self) -> Vec<i64> {
        let mut rng = rand::thread_rng();
        let mut omega_k: Vec<i64> = Vec::new();
        for i in 0..256 {
            omega_k.push(rng.gen_range(0..q));
        }
        self.omega_k = omega_k;
        omega_k
    }


    pub fn verify_b_prime_prime(&self, b_prime_prime_k : Polynomial<i64>) -> bool {
        
        let mut sum = 0;
        for i in 0..L {
            sum += self.phi_k[i] * b_prime;
        }

        // check that the constant term is equal to the above stuff.
        poly.eval(0) == check_val
    }


    pub fn sample_jl_projection() -> Array2<i64> {

        let between = Uniform::from(-1..=1);

        let mut rng = rand::thread_rng();

        let mut pi_i : Array2<i64> = Array2::zeros((256, N));

        for ((i, j), value) in pi_i.indexed_iter_mut() {
            *value = between.sample(&mut rng);
            //println!("matrix[{}][{}] = {}", i, j, value);
        }
        // TODO store pi_i in the verifier's data
        Pi.push(pi_i);

        pi_i
    }


    pub fn valid_projection(projection: Array2<Polynomial<i64>>) -> bool {
        self.projection = projection;
        let val : f64 = 128.;
        let total_norm = compute_total_norm(projection);
        println!("TOTAL NORM OF JL PROJECTION: {}, {}",total_norm, val.sqrt()*(BETA_BOUND as f64));
        return total_norm <= val.sqrt()*(BETA_BOUND as f64);
    }

}

