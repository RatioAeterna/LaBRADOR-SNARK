use ndarray::{Array2, Ix2, concatenate, Axis};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;
use crate::verification::*;
use crate::structs::*;


pub struct Prover {

    // prover fields
    S: Array2<Polynomial<i64>>,

}

impl Prover {


    pub fn new(S: Array2<Polynomial<i64>>) -> Self {
        Prover {
            S,
        }
    }

    pub fn proof_gen(&self, verifier : Verifier) {
        // random generation of the polynomial matrix A
        let mut A = Array2::from_elem((R,R), Polynomial::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let a_ij = generate_polynomial(Q,D);
                if A[[i,j]] == Polynomial::new(vec![]) {
                    A[[i,j]] = a_ij.clone();
                    A[[j,i]] = a_ij;
                }
            }
        }
        println!("Generated A!");
        for row in A.rows() {
            for poly in row {
                println!("{}", poly.pretty("x"));
            }
        }

        // random generation of random polynomial matrix Phi
        let mut Phi = Array2::from_elem((N,R), Polynomial::new(vec![])); 
        for i in 0..N {
            for j in 0..R {
                let phi_ij = generate_polynomial(Q,D);
                Phi[[i,j]] = phi_ij;
            }
        }

        // Next, we need to compute 'b' such that this relation is equal to zero
        // So we compute the values of the relation and work backwards
        let mut a_product : Polynomial<i64> = Polynomial::new(vec![]);
        for i in 1..R {
            for j in 1..R {
                let inner_prod = polynomial_vec_inner_product(self.S.column(i).to_vec(), self.S.column(j).to_vec());
                let prod = A[[i,j]].clone() * inner_prod;
                a_product = a_product.clone() + prod;
            }
        }

        let mut phi_product : Polynomial<i64> = Polynomial::new(vec![]);
        for i in 1..R {
            let inner_prod = polynomial_vec_inner_product(Phi.column(i).to_vec(), self.S.column(i).to_vec());
            phi_product = phi_product.clone() + inner_prod;
        }

        let b : Polynomial<i64> = a_product + phi_product;
        println!("Generated b!");
        println!("{}", b.pretty("x"));


        //let t = Array2::<i64>::zeros((M,R)); // Inner commitment vector 
        //let t = Array2<Polynomial<i64>>::zeros((R,N)); // Inner commitment vector ( After looking at this some more, I think this is the actual dimension. Unclear, however. )
        let zero_poly = Polynomial::new(vec![0i64]);
        let mut t = Array2::from_elem((R,N), zero_poly);

        // Compute inner Ajtai commitments
        // t_i = As_i \in R_q^\kappa (should just be m-dimensional)
        for i in 0..N {
            //let t_i = A.clone().dot(S.t().column(i)); // we need to compute the transpose of S for this product to work (for some reason) // we need to compute the transpose of S for this product to work (for some reason)
            let t_i = polynomial_matrix_product(A.clone(), self.S.t().column(i).to_owned().insert_axis(Axis(1))); // we need to compute the transpose of S for this product to work (for some reason) // we need to compute the transpose of S for this product to work (for some reason)
            println!("A dim {:?}", A.dim());
            println!("S.t() dim {:?}", self.S.t().dim());
            println!("S col dim {:?}", self.S.t().column(i).dim());
            println!("t_i dim: {:?}", t_i.dim());
            t.column_mut(i).assign(&t_i.remove_axis(Axis(1)));
        }
        println!("Computed Inner Ajtai Commitments!");

        // Next, compute JL projection
        while !verifier.valid_projection(jl_project(self.S.clone())) {
            println!("Verifier rejected projection!");
        }
        println!("JL Projection complete and accepted");

        // First Aggregation step
        

        // TODO the paper mentions 'log' but it's not AT ALL clear whether this is base 2, 10, e, etc.
        // We will go with 10 for now.

        let upper_bound : usize = (128.0f64 / (q as f64).log10()).ceil();
        for k in 1..upper_bound {
            let phi_k = verifier.generate_phi_k();
            let omega_k = verifier.generate_omega_k();


            // NOTE: for our first basic rudimentary implementation, consider that
            // a_prime_ij is just going to be a_ij... given that our F' is just going to be F since
            // they both have zero constant coefficient. Simple F' for now.
            for i in 0..r {
                for j in 0..r {
                    let a_prime_ij : Polynomial<i64> = A[[i,j]];
                    let a_prime_prime_ij = multiply_poly_ints(a_prime_ij, phi_k);

                    // TODO is the Phi[i] indexing the right way in terms of column vs row? Just
                    // .col()
                    let phi_prime_i : Vec<Polynomial<i64>> = Phi[i].to_vec();

                    // side note: consider that pi_i^(j) is the jth row of Pi_i for j = 1, ..., 256.

                    let Pi_i = verifier.get_Pi_i();

                    let mut rhs : Polynomial<i64> = Polynomial::new(vec![]);

                    // TODO I think this should be 0..256 or something depending on how for loops
                    // work
                    for j in 1..256 {
                        let conj = sigma_inv(Pi_i.row(j));
                        let omega_k_j = omega_k[j];

                        let res = scale_polynomial(conj, omega_k_j);
                        rhs += res;
                    }

                    let phi_i_prime_prime = multiply_poly_vec_ints(phi_prime_i , phi_k) + rhs; 
                }
            }

            verifier.verify_b_prime_prime(b_prime_prime_k);
        }



        // TODO fill the proof transcript with all the relevant data and return
        let proof_transcript : Transcript = Transcript { };
        proof_transcript
    }

    pub fn jl_project(&self) -> Array2<Polynomial<i64>> {
        // verifier sends random matrices in {-1,0,1}
        // TODO/NOTE: there is great confusion in the paper.
        // they state that each projection matrix pi_i generated by the verifier is of dimension
        // 256 x (n*d).. which doesn't really make sense because s_i is of dimension n x 1. Strange.
        // We'll ignore the factor of d for now and potentially come back to it.
        let mut projection : Array2<Polynomial<i64>> = Array2::zeros((256,1));
        for i in 0..R {
            let pi_i = sample_jl_projection();

            let mut product : Array2<Polynomial<i64>> = Array2::zeros((256,1));
            for row in 0..pi_i.nrows() {
                let mut sum : Polynomial<i64> = Polynomial::new(vec![]);
                for col in 0..pi_i.ncols() {
                    // scale the polynomial by either -1, 0, or 1
                    sum = sum + scale_polynomial(self.S.column(i)[col].clone(), pi_i[[row,col]] as f32);
                }
                product[[row, 0]] = sum;
            }
            projection = projection + product;
        }
        projection
    }

}

// generates and returns the SIS vector s
pub fn generate_witness() -> Array2<Polynomial<i64>> {
    
    let mut S = Array2::from_elem(Ix2(N, R), Polynomial::new(vec![])); 

    // Random Generation of Witness matrix S
    
    // total summed norm:
    let mut norm_sum: i64 = 0;


    for i in 0..N {
        for j in 0..R {
            let mut s_ij = generate_polynomial(Q, D);
            S[[i,j]] = s_ij.clone();
            // add norm to norm sum
            norm_sum += poly_norm(s_ij) as i64;
        }
    }

    // if too big, scale down:
    if norm_sum > (i64::pow(BETA_BOUND,2)) {
        let scale_factor : f32 = ((i64::pow(BETA_BOUND,2)) / norm_sum) as f32;
        // scale each polynomial in the matrix by scale factor
        for i in 0..N {
            for j in 0..R {
                S[[i,j]] = scale_polynomial(S[[i,j]].clone(), scale_factor);
            }
        }
    }

    S
}
