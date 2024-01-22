use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;
use crate::structs::*;

pub struct Verifier {

    //projection : Option<&'a Array2<Polynomial<i64>>>,
    //psi_k : Option<Vec<i64>>,
    //omega_k : Option<&'a Vec<i64>>,
    b_prime : Option<i64>,

    // JL projection matrix
    Pi : Option<Vec<Array2<i64>>>,
}

impl Verifier {

    pub fn new() -> Self {
        Verifier {
            //projection: None,
            //psi_k: None,
            //omega_k: None,
            b_prime: None,
            Pi: None,
        }
    }


    pub fn verify(&self, st: &State, proof : Transcript, crs : &CRS) -> bool {
        
        // CHECK 8
        // check that g_{ij} == g_{ji} i.e., matrix Gij is symmetric
        // TODO is it faster to only check some values? I really don't think this makes a
        // difference.
        for i in 0..R {
            for j in 0..R {
                if (proof.Gij[[i,j]] != proof.Gij[[j,i]]) {
                    return false;
                }
            }
        }



        // CHECK 16
        let lhs : Polynomial<i64> = polynomial_vec_inner_product(&proof.z, &proof.z);
        let mut rhs : Polynomial<i64> = Polynomial::new(vec![]);

        for i in 0..R {
            for j in 0..R {
                rhs = rhs + (&proof.Gij[[i,j]] * &proof.c[i] * &proof.c[j]);
            }
        }
        if (lhs != rhs) {
            return false;
        }

        true
    }

    
    pub fn get_Pi_i(&self, i : usize) -> &Array2<i64> {
        &self.Pi.as_ref().unwrap()[i]
    }


    // TODO do we want to STORE alpha, beta in the Verifier struct?
    pub fn fetch_alpha(&self) -> Vec<Polynomial<i64>> {
        let mut alpha = vec![]; 
        for i in 0..K {
            alpha.push(generate_polynomial(Q,D));
        }
        alpha
    }

    pub fn fetch_beta(&self) -> Vec<Polynomial<i64>> {
        let mut beta = vec![]; 
        let upper_bound : usize = (128.0f64 / (Q as f64).log10()).ceil() as usize;
        for i in 0..upper_bound {
            beta.push(generate_polynomial(Q,D));
        }
        beta
    }


    // fetch a challenge polynomial from the challenge space \mathcal{C} satisfying a number of
    // criteria
    pub fn fetch_challenge(&self) -> Polynomial<i64> {
        // particular challenge coefficient distribution described on page 6.
        let mut coeff_dist : Vec<i64> = vec![];
        for i in 0..23 {
            coeff_dist.push(0);
        }
        for i in 0..31 {
            coeff_dist.push(1);
        }
        for i in 0..10 {
            coeff_dist.push(2);
        }

        let candidate = generate_polynomial_picky(Q,D as usize, coeff_dist.clone());
        // TODO... I don't think this norm should be squared, as it is.. which would perhaps give you 71^2.. definitely fix this if
        // needed.
        assert!(poly_norm(&candidate) == TAU, "Incorrect l2 norm of challenge polynomial");
    
        while operator_norm(&candidate) > T {
            assert!(poly_norm(&candidate) == TAU, "Incorrect l2 norm of challenge polynomial");
            let candidate = generate_polynomial_picky(Q,D as usize, coeff_dist.clone());
        }
        candidate
    }

    pub fn generate_psi(&self) -> Vec<i64> {
        let mut rng = rand::thread_rng();
        let mut psi_k: Vec<i64> = Vec::new();
        for i in 0..L {
            psi_k.push(rng.gen_range(0..Q));
        }
        psi_k
    }
            

    pub fn generate_omega(&self) -> Vec<i64> {
        let mut rng = rand::thread_rng();
        let mut omega_k: Vec<i64> = Vec::new();
        for i in 0..256 {
            omega_k.push(rng.gen_range(0..Q));
        }
        //self.omega_k = Some(&omega_k);
        omega_k 
    }


    /*
    pub fn verify_b_prime_prime(&self, b_prime_prime_k : Polynomial<i64>) -> bool {
        // TODO again column vs row not sure.
        // Also self vs no self keyword not sure.
        let prod = vec_inner_product(self.omega_k.unwrap(), self.projection.unwrap().column(0).to_vec());
        let mut sum = 0;
        for i in 0..L {
            sum += self.psi_k.unwrap()[i] * self.b_prime.unwrap();
        }
        let check_val = prod + sum;

        // check that the constant term is equal to the above stuff.
        b_prime_prime_k.eval(0) == check_val
    }
    */


    pub fn sample_jl_projection(&mut self) -> &Array2<i64> {

        let between = Uniform::from(-1..=1);

        let mut rng = rand::thread_rng();

        let mut Pi_i : Array2<i64> = Array2::zeros((256, N*(D as usize)));

        for ((i, j), value) in Pi_i.indexed_iter_mut() {
            *value = between.sample(&mut rng);
            //println!("matrix[{}][{}] = {}", i, j, value);
        }
        // TODO store pi_i in the verifier's data
        self.Pi.as_mut().unwrap().push(Pi_i);

        //&(self.Pi.unwrap().last().unwrap())
        // TODO don't entirely understand the functionality of this line.. but seems to work.
        let Pi_ref = self.Pi.as_ref().and_then(|pi| pi.last()).unwrap();
        Pi_ref
    }


    pub fn valid_projection(&mut self, projection: &Vec<i64>) -> bool {
        let val : f64 = 128.;
        let norm = l2_norm(projection);
        println!("TOTAL NORM OF JL PROJECTION: {}, {}", norm, val.sqrt()*(BETA_BOUND as f64));
        return norm <= val.sqrt()*(BETA_BOUND as f64);
    }

}

