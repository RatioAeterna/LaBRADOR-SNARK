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


    // TODO do we want to STORE alpha, beta in the Verifier struct?
    pub fn fetch_alpha(&self) -> Vec<Polynomial<i64>> {
        let mut alpha = vec![]; 
        for i in 0..K {
            alpha.push(generate_polynomial(q,d));
        }
        alpha
    }

    pub fn fetch_beta(&self) -> Vec<Polynomial<i64>> {
        let mut beta = vec![]; 
        let upper_bound : usize = (128.0f64 / (q as f64).log10()).ceil();
        for i in 0..upper_bound {
            beta.push(generate_polynomial(q,d));
        }
        beta
    }


    // fetch a challenge polynomial from the challenge space \mathcal{C} satisfying a number of
    // criteria
    pub fn fetch_challenge(&self) -> Polynomial<i64> {
        // particular challenge coefficient distribution described on page 6.
        let coeff_dist : Vec<i64> = vec![];
        for i in 0..23 {
            coeff_dist.push(0);
        }
        for i in 0..31 {
            coeff_dist.push(1);
        }
        for i in 0..10 {
            coeff_dist.push(2);
        }

        let candidate = generate_polynomial_picky(q,d, coeff_dist);
        // TODO... I don't think this norm should be squared, as it is.. which would perhaps give you 71^2.. definitely fix this if
        // needed.
        assert!(poly_norm(candidate) == 71, "Incorrect l2 norm of challenge polynomial");
    
        while operator_norm(candidate) > 15 {
            assert!(poly_norm(candidate) == 71, "Incorrect l2 norm of challenge polynomial");
            let candidate = generate_polynomial_picky(q,d, coeff_dist);
        }
        candidate
    }

    pub fn generate_phi(&self) -> Vec<i64> {
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
        // TODO again column vs row not sure.
        // Also self vs no self keyword not sure.
        let prod = vec_inner_product(self.omega_k, self.projection.column(0).to_vec());
        let mut sum = 0;
        for i in 0..L {
            sum += self.phi_k[i] * b_prime;
        }
        let check_val = prod + sum;

        // check that the constant term is equal to the above stuff.
        b_prime_prime_k.eval(0) == check_val
    }


    pub fn sample_jl_projection(&self) -> Array2<i64> {

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


    pub fn valid_projection(&self, projection: Array2<Polynomial<i64>>) -> bool {
        self.projection = projection;
        let val : f64 = 128.;
        let total_norm = compute_total_norm(projection);
        println!("TOTAL NORM OF JL PROJECTION: {}, {}",total_norm, val.sqrt()*(BETA_BOUND as f64));
        return total_norm <= val.sqrt()*(BETA_BOUND as f64);
    }

}

