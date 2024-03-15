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
use std::sync::atomic::{AtomicBool, Ordering};
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive};


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
        }
        println!("Computed Inner Ajtai Commitments!");


        println!("Generating Gij...");
        // NOTE: This is a symmetric matrix
        let mut Gij = Array2::from_elem((R,R), R_q::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let col_i = self.S.column(i).to_vec();
                let col_j = self.S.column(j).to_vec();

                let res : R_q = polynomial_vec_inner_product(&col_i, &col_j);
                if Gij[[i,j]] == R_q::new(vec![]) {
                    Gij[[i,j]] = res
                }
            }
        }
        println!("Finished generating Gij.");

        println!("Generating polynomial decomposition of t_i");
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

        println!("Generating polynomial decomposition of g_i");
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

        println!("Generating u_1...");
        // TODO investigate weird discrepancy between KAPPA_1 and KAPPA_2 dimension... 
        // FOR NOW, we assume they are equal rank. But I *think* you zero pad rhs if this isn't the
        // case.
        let u_1 : Vec<R_q> = add_poly_vec(&lhs, &rhs);
        //println!("IMPORTANT PRINT!!! lhs: {:?}, rhs: {:?}", &lhs, &rhs);

        println!("Computing JL projection...");
        // Next, compute JL projection
        let mut projection_ints = self.jl_project();
        println!("Computed projection..");
        let mut rejections = 0;
        while !self.verifier.valid_projection(&projection_ints) {
            println!("Verifier rejected projection!");
            rejections += 1;
            if(rejections > 2) {panic!("failed JL...");}
            projection_ints = self.jl_project();
        }
        println!("JL Projection complete and accepted");

        let mut projection = Z_q::lift(&projection_ints);

        // First Aggregation step

        // TODO the paper mentions 'log' but it's not AT ALL clear whether this is base 2, 10, e, etc.
        // We will go with 10 for now.
        let upper_bound : usize = std::cmp::min(K, (128.0f64 / (*Q as f64).log2()).ceil() as usize);


        let mut psi : Vec<Vec<Z_q>> = vec![];
        let mut omega : Vec<Vec<Z_q>> = vec![];


        // TODO these comments aren't precise enough in dimensions
        // vector containing all phi'' for k in 1..(upper_bound), i in 0..(R-1)
        // TODO also perhaps switch this to an actual Vec<Array2<R_q>> later.
        let mut phi_prime_prime_vec : Vec<Vec<Vec<R_q>>> = vec![];

        // vector containing all b'' for k in 1..(upper_bound)
        let mut b_prime_prime_vec : Vec<R_q> = vec![];


        println!("First aggregation step...");
        for k in 1..(upper_bound+1) {
            println!("k={}, upper_bound: {}", k, upper_bound);
            let psi_k = self.verifier.generate_psi();
            let omega_k = self.verifier.generate_omega();

            let mut phi_prime_prime_k : Vec<Vec<R_q>> = vec![];

            psi.push(psi_k.clone());
            omega.push(omega_k);

            // NOTE: for our first basic rudimentary implementation, consider that
            // a_prime_ij is just going to be a_ij... given that our F' is just going to be F since
            // they both have zero constant coefficient. Simple F' for now.

            // TODO obviously we'll iterate the below look across all L in the end,
            // AND a_prime_prime will have a (k) superscript... so this will be a full K length vec eventually
            let mut a_prime_prime : Array2<R_q> = Array2::zeros((R,R));

            for i in 0..R {
                let mut phi_i_prime_prime : Vec<R_q> = vec![];
                for j in 0..R {

                    let a_prime_ij : &R_q = &st.a_prime_k[0][[i,j]];
                    let a_prime_prime_ij : R_q = multiply_poly_ints(&a_prime_ij, &psi[k-1]);
                    a_prime_prime[[i,j]] = a_prime_prime_ij;
                }

                // TODO is the phi_k[0][i] indexing the right way in terms of column vs row? Just
                // .col()
                // TODO we're taking the k=0 index for now.. disregarding the rest
                let phi_prime_i : Vec<R_q> = st.phi_k[0].column(i).to_vec();

                let mut lhs : Vec<R_q> = gen_empty_poly_vec(N);
                for l in 0..L {                                 //NOTE: This is "psi_k"
                    lhs = add_poly_vec(&multiply_poly_vec_ints(&phi_prime_i , &psi.last().unwrap()), &lhs);
                }

                // side note: consider that pi_i^(j) is the jth row of Pi_i for j = 1, ..., 256.
                let Pi_i = self.verifier.get_Pi_i(i);
                let mut rhs : Vec<R_q> = vec![R_q::new(vec![]); N as usize];
                for j in 0..256 {
                    let bolded_pi_poly_vec : Vec<R_q> = concat_coeff_reduction(&Z_q::lift(&Pi_i.row(j).to_vec()));
                    //let bolded_pi_poly_vec : Vec<R_q> = concat_coeff_reduction(&Pi_i.row(j).to_vec());
                    let conj = sigma_inv_vec(&bolded_pi_poly_vec);
                    let omega_k_j = &omega[0][j];
                    let res = scale_poly_vec(&conj, f32::from(omega_k_j));
                    rhs = add_poly_vec(&rhs, &res);
                }
                phi_i_prime_prime = add_poly_vec(&lhs, &rhs); 
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
            let b_prime_prime_k = lhs.clone() + rhs.clone();

            // FIRST, we sanity check to see if this whole new "Aggregation" f''^(k) relation is zero..
            // TODO debugging, not part of the protocol

            let mut b_prime_prime_0 : Z_q = self.verifier.fetch_alleged_b_prime_prime_cc(&omega.last().unwrap(), &psi.last().unwrap(), &projection);


            // TODO DEBUGGING: Check left side of equation
            // ===========================================
            let mut b_prime_prime_test : R_q = R_q::zero();
            let mut b_prime : R_q = R_q::zero();

            for i in 0..R {
                for j in 0..R {
                    let a_prime_ij : &R_q = &st.a_prime_k[0][[i,j]];
                    let prod = polynomial_vec_inner_product(&self.S.column(i).to_vec(), &self.S.column(j).to_vec());
                    let res = a_prime_ij * &prod;
                    b_prime = b_prime + res;
                }
                let phi_prime_i : Vec<R_q> = st.phi_k[0].column(i).to_vec();
                b_prime = b_prime + polynomial_vec_inner_product(&phi_prime_i, &self.S.column(i).to_vec());
            }
            println!("b_prime... {}", &b_prime);

            b_prime_prime_test = multiply_poly_ints(&b_prime, &psi_k);


            // recompute by multiplying by psi..
            let mut b_prime_prime_candidate : R_q = R_q::zero();
            for i in 0..R {
                for j in 0..R {
                    let a_prime_ij : &R_q = &st.a_prime_k[0][[i,j]];
                    let prod = polynomial_vec_inner_product(&self.S.column(i).to_vec(), &self.S.column(j).to_vec());
                    let res = a_prime_ij * &prod;
                    let resres = multiply_poly_ints(&res, &psi_k);
                    b_prime_prime_candidate = b_prime_prime_candidate + res;
                }
                let phi_prime_i : Vec<R_q> = st.phi_k[0].column(i).to_vec();
                let prod2 = polynomial_vec_inner_product(&phi_prime_i, &self.S.column(i).to_vec());
                b_prime_prime_candidate = b_prime_prime_candidate + multiply_poly_ints(&prod2, &psi_k);
            }
            assert!(b_prime_prime_test == b_prime_prime_candidate, "ASSERT FAILED... multiply by psi after: {}, multiply during: {}", b_prime_prime_test, b_prime_prime_candidate);




            // TODO DEBUG: CHECK CONJUGATION AUTOMORPHISM INVARIANT 
            /*
            for j in 0..256 {
                //let mut testres : R_q = R_q::zero();
                let mut p_j : Z_q = Z_q::from(0);
                let mut right_side : R_q = R_q::zero();

                for i in 0..R {
                    // compute left side of invariant..
                    println!("\n\ni = {}", i);
                    let Pi_i = self.verifier.get_Pi_i(i);
                    let s_i_coeffs : Vec<Z_q> = witness_coeff_concat(&self.S.column(i).to_vec());
                    println!("printing s_i coeffs: {:?}... and here's actual s_i!: {:?}", s_i_coeffs, &self.S.column(i).to_vec());
                    println!("printing pi_i^(j) coeffs: {:?}...", &Pi_i.row(j).to_vec());

                    p_j = p_j + vec_inner_product_Z_q(&Pi_i.row(j).to_vec(), &s_i_coeffs);

                   
                    //let bolded_pi_poly_vec : Vec<R_q> = concat_coeff_reduction(&Z_q::lift(&Pi_i.row(j).to_vec()));
                    let bolded_pi_poly_vec : Vec<R_q> = concat_coeff_reduction(&Pi_i.row(j).to_vec());
                    let mut conj = sigma_inv_vec(&bolded_pi_poly_vec);
                    //conj = bolded_pi_poly_vec;
                    //assert!(bolded_pi_poly_vec == conj, "sigma inv failed! {:?} {:?}", bolded_pi_poly_vec, conj);
                    //conj = scale_poly_vec_int(&conj, &omega.last().unwrap()[j]);


                    let prod = polynomial_vec_inner_product(&conj, &self.S.column(i).to_vec());

                    right_side = right_side + prod;
                    assert!(right_side.eval(Z_q::from(0)) == p_j, "INVARIANT BROKEN! {} {}", right_side.eval(Z_q::from(0)), &p_j); 
                }
                assert!(p_j == projection[j], "super super broken");
                //assert!(right_side.eval(Z_q::from(0)) == p_j, "INVARIANT BROKEN! {} {}", right_side.eval(Z_q::from(0)), &p_j); 
            } 

            //let right_side_prod = Z_q::from(vec_inner_product(&Z_q::lift_inv(&omega.last().unwrap()), &projection));
            //assert!(testres.eval(Z_q::from(0)) == right_side_prod, "not equal! {} {}", testres.eval(Z_q::from(0)), right_side_prod);

            

            let real_cc = (&lhs + &rhs - R_q::new(vec![b_prime_prime_0])).eval(Z_q::from(0));
            println!("alleged b prime prime 0... {}", b_prime_prime_0);
            println!("Real constant coefficient diff: {} ", real_cc);
            */



            self.verifier.verify_b_prime_prime(&b_prime_prime_k, &omega.last().unwrap(), &psi.last().unwrap(), &projection);

            phi_prime_prime_vec.push(phi_prime_prime_k);
            b_prime_prime_vec.push(b_prime_prime_k);
        }

        println!("fetching alpha and beta...");
        
        let alpha : Vec<R_q> = self.verifier.fetch_alpha();
        let beta : Vec<R_q> = self.verifier.fetch_beta();

        println!("Generating phi final..");
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

        println!("Generating Hij..");
        // NOTE: This is *also* a symmetric matrix
        let mut Hij = Array2::from_elem((R,R), R_q::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                //println!("i : {}, j: {}", i, j);
                let s_i = self.S.column(i).to_vec();
                let s_j = self.S.column(j).to_vec();
                //println!("s_i: {:?}, s_j: {:?}", s_i, s_j);

                let phi_i = &phi_final[i];
                let phi_j = &phi_final[j];
                //println!("phi_i: {:?}, phi_j: {:?}", phi_i, phi_j);

                MOD_SUSPENSION.store(true, Ordering::SeqCst);

                let sum = polynomial_vec_inner_product(phi_i, &s_j) + polynomial_vec_inner_product(phi_j, &s_i);

                //println!("sum: {}", sum);

                // NOTE here we're computing the 
                // multiplicative inverse of 2, i.e., "1/2"
                let two = BigInt::from(2);
                let q_bigint = BigInt::from(*Q);
                let inv : i128 = two.modpow(&(q_bigint.clone() - BigInt::from(2)), &q_bigint).to_i128().unwrap();


                let mut res = scale_polynomial_int(&sum, &Z_q::from(inv));

                //let mut res = scale_polynomial_rational(&sum, &Z_q::from(1), &Z_q::from(2));
                //println!("sum/2: {}", res);

                MOD_SUSPENSION.store(false, Ordering::SeqCst);
                res = res.recompute_mod();

                //println!("res: {}", res);

                if Hij[[i,j]] == R_q::new(vec![]) {
                    Hij[[i,j]] = res;
                }
            }
        }

        /*
        println!("Printing Hij!");
        for row in Hij.rows() {
            for poly in row {
                print!("{}, ", poly.eval(Z_q::new(0)));
            }
            println!("\n");
        }
        */




        println!("Computing u_2...");
        // Compute u_2
        let mut u_2 : Vec<R_q> = vec![R_q::new(vec![]); KAPPA_2 as usize];
        for i in 0..R {
            for j in i..R {
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

        let mut z : Vec<R_q> = vec![R_q::new(vec![]); N as usize];
        let mut c_vec : Vec<R_q> = vec![];

        println!("challenge polynomial fetching..");
        //println!("JUST FOR SANITY, WITNESS SHAPE: {:?}", &self.S.shape()); 
        for i in 0..R {
            let c_i = self.verifier.fetch_challenge();
            c_vec.push(c_i);
            let prod : Vec<R_q> = poly_by_poly_vec(c_vec.last().as_ref().unwrap(), &self.S.column(i).to_vec());
            //println!("WITNESS COL VEC SIZE: {}", &self.S.column(i).to_vec().len()); 
            //println!("PROD VEC SIZE: {}", &prod.len()); 
            for j in 0..N {
                z[j] = &z[j] + &prod[j];
            }
        }
        //println!("z length.. {}" , &z.len()); 
        println!("Filling proof transcript.");
        // TODO fill the proof transcript with all the relevant data and return
        let proof_transcript : Transcript = Transcript { u_1, projection, psi, omega, b_prime_prime : b_prime_prime_vec, alpha, beta, u_2, c : c_vec, z, t_i_all, Gij, Hij };
        proof_transcript
    }

    pub fn jl_project(&mut self) -> Vec<i128> {
        // verifier sends random matrices in {-1,0,1}... 
        let mut projection : Vec<i128> = vec![0 ; 256];
        for i in 0..R {
            let Pi_i : &Array2<i128> = self.verifier.sample_jl_projection();
            //println!("Got Pi_i for i={}",i);
            let s_i_coeffs : Array2<i128> = vec_to_column_array(&Z_q::lift_inv(&witness_coeff_concat(&self.S.column(i).to_vec())));
            //let s_i_coeffs : Array2<Z_q> = vec_to_column_array(&witness_coeff_concat(&self.S.column(i).to_vec()));
            //println!("Got s_i coeffs");
            // NOTE: for reference, this is a 256x(ND) multiplied by an (ND)x1, giving a 256x1
            // which we turn into a vec
            let product = matmul(Pi_i, &s_i_coeffs).column(0).to_vec();
            //println!("computed product");
            projection = add_vecs(&projection, &product);
        }
        projection
    }
}

// generates and returns the SIS vector s
pub fn generate_witness() -> Array2<R_q> {
    
    let mut S = Array2::from_elem(Ix2(N, R), R_q::new(vec![])); 
    let mut rng = rand::thread_rng();

    // Random Generation of Witness matrix S
    
    // total summed norm:
    let mut norm_sum: i128 = 0;

    for i in 0..N {
        for j in 0..R {
            let mut s_ij = generate_polynomial(*Q, D);
            // TODO fix clone
            S[[i,j]] = s_ij.clone();
            // add norm to norm sum
            norm_sum += poly_norm(&s_ij) as i128;
        }
    }

    /*
    for i in 0..N {
        for j in 0..R {
            println!("val at i={}, j={}: {}", i, j, S[[i,j]]);
        }
    }
    */



    // if too big, scale down:
    while norm_sum > (i128::pow(*BETA_BOUND,2)) {
        println!("norm sum too big! {} {}", norm_sum, (*BETA_BOUND)*(*BETA_BOUND));

        let i = rng.gen_range(0..N);
        let j = rng.gen_range(0..R);

        let old_poly_norm = poly_norm(&S[[i,j]]) as i128;

        //println!("i {}, j {},.... {}", i, j, poly_norm(&S[[i,j]]) as i128);
        //println!("POLY NORM: {}", poly_norm(&S[[i,j]]) as i128);
        //println!("norm sum: {}", norm_sum);
        //println!("Old Sij: {}", &S[[i,j]]);

        let mut s_ij = reduce_polynomial(&S[[i,j]]);

        //println!("New Sij: {}", &s_ij);
        let new_poly_norm = poly_norm(&s_ij) as i128;
        //println!("norm diff... {}", old_poly_norm - new_poly_norm);
        //let mut s_ij = R_q::zero();

        S[[i,j]] = s_ij.clone();

        norm_sum -= (old_poly_norm - new_poly_norm);
    }
    //println!("SUPER SNITY CHECK VAL!! {}", S[[0,0]].clone());
    println!("norm of witness! squared version: {}, sqrt version: {}", norm_sum, (norm_sum as f64).sqrt() as i128);
    S
}
