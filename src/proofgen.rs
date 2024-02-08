use ndarray::{Array2, Ix2, concatenate, Axis}; 
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::algebraic::*;
use crate::util::*;
use crate::constants::*;
use crate::verification::*;
use crate::structs::*;
use num_traits::Zero;


pub struct Prover<'a> {

    // prover fields
    S: &'a Array2<R_q>,

    // reference to verifier
    // TODO does this even need to be mutable????
    verifier: &'a mut Verifier,

}

impl<'a> Prover<'a> {


    pub fn new(S: &'a Array2<R_q>, verifier: &'a mut Verifier) -> Self {
        Prover {
            S,
            verifier,
        }
    }

    pub fn proof_gen(&mut self, st : &State, crs : &CRS) -> Transcript {
        let mut t_i_all : Vec<Vec<R_q>> = vec![];

        println!("Generating inner Ajtai commitments...");
        // Compute inner Ajtai commitments
        // t_i = As_i \in R_q^\kappa (should just be m-dimensional)
        for i in 0..R {
            // TODO this doesn't look nice, fix later
            let column = self.S.column(i).to_owned(); // Convert column view to an owned 1D array
            // TODO this is also bad because you're hardcoding N for now. Fix later.
            let column_as_2d = column.into_shape((N as usize, 1)).unwrap(); // Reshape into 2D array with one column
            let t_i = polynomial_matrix_product(&crs.A_mat, &column_as_2d).column(0).to_vec(); 
            t_i_all.push(t_i);
            /*
            println!("A dim {:?}", crs.A_mat.dim());
            println!("S.t() dim {:?}", self.S.t().dim());
            println!("S col dim {:?}", self.S.t().column(i).dim());
            println!("t_i dim: {:?}", t_i.dim());
            */
        }
        println!("Computed Inner Ajtai Commitments!");


        // NOTE: This is a symmetric matrix
        let mut Gij = Array2::from_elem((R,R), R_q::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let col_i = self.S.t().column(i).to_vec();
                let col_j = self.S.t().column(j).to_vec();

                let res : R_q = polynomial_vec_inner_product(&col_i, &col_j);
                if Gij[[i,j]] == R_q::new(vec![]) {
                    Gij[[i,j]] = res;
                }
            }
        }

        let mut lhs = gen_empty_poly_vec(KAPPA_1 as usize);
        for i in 0..R {
            let t_i_decomposed : Vec<Vec<R_q>> = decompose_polynomial_vec(&t_i_all[i], *B_1, *T_1);
            for k in 0..(*T_1 as usize) {
                let B_ik = crs.B_mat.get(&(i,k)).unwrap();
                let t_mat = vec_to_column_array(&t_i_decomposed[k]);
                // NOTE this line might look somewhat confusing, but consider that the matrix
                // product is just a KAPPA_1 x 1 Array2 (which we turn into a vec by taking that
                // first column
                let prod = polynomial_matrix_product(B_ik, &t_mat).column(0).to_vec();
                lhs = add_poly_vec(&prod, &lhs);
            }
        }

        let mut rhs = gen_empty_poly_vec(KAPPA_2 as usize);
        for i in 0..R {
            for j in i..R {
                let g_ij : Vec<R_q> = decompose_polynomial(&Gij[[i,j]], *B_2, *T_2);
                for k in 0..(*T_2 as usize) {
                    let C_ijk = crs.C_mat.get(&(i,j,k)).unwrap().column(0).to_vec();
                    let g_ij_k : &R_q = &g_ij[k];
                    let prod = poly_by_poly_vec(g_ij_k, &C_ijk);
                    rhs = add_poly_vec(&prod, &rhs);
                }
            }
        }

        // TODO investigate weird discrepancy between KAPPA_1 and KAPPA_2 dimension... 
        // FOR NOW, we assume they are equal rank. But I *think* you zero pad rhs if this isn't the
        // case.
        let u_1 : Vec<R_q> = add_poly_vec(&lhs, &rhs);

        // Next, compute JL projection
        let mut projection = self.jl_project();
        while !self.verifier.valid_projection(&projection) {
            println!("Verifier rejected projection!");
            projection = self.jl_project();
        }
        println!("JL Projection complete and accepted");

        // First Aggregation step

        // TODO the paper mentions 'log' but it's not AT ALL clear whether this is base 2, 10, e, etc.
        // We will go with 10 for now.
        let upper_bound : usize = (128.0f64 / (Q as f64).log10()).ceil() as usize;


        let mut psi : Vec<Vec<Z_q>> = vec![];
        let mut omega : Vec<Vec<Z_q>> = vec![];


        // TODO these comments aren't precise enough in dimensions
        // vector containing all phi'' for k in 1..(upper_bound), i in 0..(R-1)
        // TODO also perhaps switch this to an actual Vec<Array2<R_q>> later.
        let mut phi_prime_prime_vec : Vec<Vec<Vec<R_q>>> = vec![];

        // vector containing all b'' for k in 1..(upper_bound)
        let mut b_prime_prime_vec : Vec<R_q> = vec![];


        for k in 1..(upper_bound+1) {
            let psi_k = self.verifier.generate_psi();
            let omega_k = self.verifier.generate_omega();

            let mut phi_prime_prime_k : Vec<Vec<R_q>> = vec![];

            psi.push(psi_k);
            omega.push(omega_k);

            // NOTE: for our first basic rudimentary implementation, consider that
            // a_prime_ij is just going to be a_ij... given that our F' is just going to be F since
            // they both have zero constant coefficient. Simple F' for now.

            // TODO obviously we'll iterate the below look across all L in the end,
            // AND a_prime_prime will have a (k) superscript... so this will be a full K length vec eventually
            let mut a_prime_prime : Array2<R_q> = Array2::zeros((R,R));

            for i in 0..R {
                let phi_i_prime_prime : Vec<R_q> = vec![];
                for j in 0..R {
                    let a_prime_ij : &R_q = &st.a_prime_k[0][[i,j]];
                    let a_prime_prime_ij : R_q = multiply_poly_ints(&a_prime_ij, &psi.last().unwrap());
                    a_prime_prime[[i,j]] = a_prime_prime_ij;

                    // TODO is the phi_k[0][i] indexing the right way in terms of column vs row? Just
                    // .col()
                    // TODO we're taking the k=0 index for now.. disregarding the rest
                    let phi_prime_i : Vec<R_q> = st.phi_k[0].column(i).to_vec();


                    let mut lhs : Vec<R_q> = gen_empty_poly_vec(N);
                    for l in 1..(L+1) {                                 //NOTE: This is "psi_k"
                        lhs = add_poly_vec(&multiply_poly_vec_ints(&phi_prime_i , &psi.last().unwrap()), &lhs);
                    }

                    // side note: consider that pi_i^(j) is the jth row of Pi_i for j = 1, ..., 256.
                    let Pi_i = self.verifier.get_Pi_i(i);

                    let mut rhs : Vec<R_q> = vec![R_q::new(vec![]); N as usize];

                    for j in 0..256 {
                        let bolded_pi_poly_vec : Vec<R_q> = concat_coeff_reduction(&Z_q::lift(&Pi_i.row(j).to_vec()));
                        let conj = sigma_inv_vec(&bolded_pi_poly_vec);
                        let omega_k_j = omega.last().unwrap()[j];

                        //let res = scale_polynomial(&conj, omega_k_j as f32);
                        //rhs = rhs + res;
    
                        let res = scale_poly_vec(&conj, f32::from(omega_k_j));
                        rhs = add_poly_vec(&rhs, &res);
                    }

                    let phi_i_prime_prime = add_poly_vec(&lhs, &rhs); 
                }
                phi_prime_prime_k.push(phi_i_prime_prime);
            }

            let mut lhs = R_q::new(vec![Z_q::zero()]);
            for i in 0..R {
                for j in 0..R {
                    // TODO I think this is S column, but might be row. Double check later.
                    let prod = polynomial_vec_inner_product(&self.S.column(i).to_vec(), &self.S.column(j).to_vec());
                    let res = &a_prime_prime[[i,j]] * prod;
                    lhs = lhs + res;
                }
            }

            let mut rhs = R_q::new(vec![Z_q::zero()]);
            for i in 0..R {
                let res = polynomial_vec_inner_product(&phi_prime_prime_k[i], &self.S.column(i).to_vec());
                rhs = rhs + res;
            }
            let b_prime_prime_k = lhs + rhs;
            self.verifier.verify_b_prime_prime(&b_prime_prime_k, &omega.last().unwrap(), &psi.last().unwrap(), &projection);

            phi_prime_prime_vec.push(phi_prime_prime_k);
            b_prime_prime_vec.push(b_prime_prime_k);
        }

        let alpha : Vec<R_q> = self.verifier.fetch_alpha();
        let beta : Vec<R_q> = self.verifier.fetch_beta();

        let mut phi_final : Vec<Vec<R_q>> = vec![];
        for i in 0..R {
            let mut lhs : Vec<R_q> = gen_empty_poly_vec(N);
            for k in 0..K {
                let res = poly_by_poly_vec(&alpha[k], &st.phi_k[k].column(i).to_vec());
                lhs = add_poly_vec(&lhs, &res);
            }
            let mut rhs : Vec<R_q> = gen_empty_poly_vec(N);
            for k in 0..upper_bound {
                let res = poly_by_poly_vec(&beta[k], &phi_prime_prime_vec[k][i]);
                rhs = add_poly_vec(&rhs , &res);
            }
            phi_final.push(add_poly_vec(&lhs, &rhs));
        }


        // NOTE: This is *also* a symmetric matrix
        let mut Hij = Array2::from_elem((R,R), R_q::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let s_i = self.S.t().column(i).to_vec();
                let s_j = self.S.t().column(j).to_vec();

                let phi_i = &phi_final[i];
                let phi_j = &phi_final[j];

                let sum = polynomial_vec_inner_product(phi_i, &s_j) + polynomial_vec_inner_product(phi_j, &s_i);
                let res = scale_polynomial(&sum, 0.5);

                if Hij[[i,j]] == R_q::new(vec![]) {
                    Hij[[i,j]] = res;
                }
            }
        }

        // Compute u_2
        let mut u_2 : Vec<R_q> = vec![R_q::new(vec![]); KAPPA_2 as usize];
        for i in 1..(R+1) {
            for j in i..(R+1) {
                for k in 0..(*T_1 as usize) {
                    // NOTE: the column(0) looks kind of suspect when the whole thing is an Array2
                    // matrix, but it's actually a column vector, so this is just an easy way to
                    // fetch this and turn it into a vec.
                    let D_ijk_vec = crs.D_mat.get(&(i,j,k)).unwrap().column(0).to_vec();
                    // TODO I don't think we can just add as such.
                    let dec_hij = decompose_polynomial(&Hij[[i,j]], *B_1, *T_1);
                    let prod = poly_by_poly_vec(&dec_hij[k], &D_ijk_vec);
                    // NOTE prod.len() = KAPPA_2 (it should at least)
                    for idx in 0..prod.len() {
                        u_2[idx] = &u_2[idx] + &prod[idx];
                    }
                }
            }
        }

        let mut z : Vec<R_q> = vec![R_q::new(vec![]); R as usize];
        let mut c_vec : Vec<R_q> = vec![];

        for i in 0..R {
            let c_i = self.verifier.fetch_challenge();
            c_vec.push(c_i);
            let prod : Vec<R_q> = poly_by_poly_vec(c_vec.last().as_ref().unwrap(), &self.S.column(i).to_vec());
            for j in 0..R {
                z[j] = &z[j] + &prod[j];
            }
        }

        // TODO fill the proof transcript with all the relevant data and return
        let proof_transcript : Transcript = Transcript { u_1, projection, psi, omega, b_prime_prime : b_prime_prime_vec, alpha, beta, u_2, c : c_vec, z, t_i_all, Gij, Hij };
        proof_transcript
    }

    pub fn jl_project(&mut self) -> Vec<i128> {
        // verifier sends random matrices in {-1,0,1}
        let mut projection : Vec<i128> = vec![0 ; 256];
        for i in 0..R {
            let Pi_i = self.verifier.sample_jl_projection();
            let s_i_coeffs : Array2<i128> = vec_to_column_array(&Z_q::lift_inv(&witness_coeff_concat(&self.S.column(i).to_vec())));
            // NOTE: for reference, this is a 256x(ND) multiplied by an (ND)x1, giving a 256x1
            // which we turn into a vec
            let product = matmul(Pi_i, &s_i_coeffs).column(0).to_vec();
            projection = add_vecs(&projection, &product);
        }
        projection
    }
}

// generates and returns the SIS vector s
pub fn generate_witness() -> Array2<R_q> {
    
    let mut S = Array2::from_elem(Ix2(N, R), R_q::new(vec![])); 

    // Random Generation of Witness matrix S
    
    // total summed norm:
    let mut norm_sum: Z_q = Z_q::zero();

    for i in 0..N {
        for j in 0..R {
            let mut s_ij = generate_sparse_polynomial(Q, D);
            // TODO fix clone
            S[[i,j]] = s_ij.clone();
            // add norm to norm sum
            norm_sum += Z_q::from(poly_norm(&s_ij));
        }
    }

    // if too big, scale down:
    if norm_sum > (i128::pow(BETA_BOUND,2)) {
        let scale_factor: f32 = (i128::pow(BETA_BOUND, 2) as f32) / (f32::from(norm_sum));
        //println!("scale factor! {}", scale_factor);
        // scale each polynomial in the matrix by scale factor
        for i in 0..N {
            for j in 0..R {
                S[[i,j]] = scale_polynomial(&S[[i,j]], scale_factor);
            }
        }
    }
    //println!("SUPER SNITY CHECK VAL!! {}", S[[0,0]].clone());
    S
}
