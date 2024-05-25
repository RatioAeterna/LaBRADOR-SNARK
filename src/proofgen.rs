use crate::algebraic::*;
use crate::constants::*;
use crate::structs::*;
use crate::util::*;
use crate::verification::*;
use ndarray::{Array2, ArrayView, Ix2};
use num_bigint::BigInt;
use num_traits::ToPrimitive;
use num_traits::Zero;
use rand::Rng;
use std::sync::atomic::Ordering;

pub struct Prover<'a> {
    // prover fields
    witness: &'a Array2<Rq>,

    // reference to verifier
    // TODO does this even need to be mutable????
    verifier: &'a Verifier<'a>,

    constants: &'a RuntimeConstants,
}

impl<'a> Prover<'a> {
    pub fn new(witness: &'a Array2<Rq>, verifier: &'a Verifier<'a>, constants: &'a RuntimeConstants) -> Self {
        Prover { witness, verifier, constants }
    }

    pub fn proof_gen(&mut self, st: &State, crs: &mut CRS) -> Transcript {
        // first, fetch witness dimension (n, r) to use later..



        let mut t_i_all: Vec<Vec<Rq>> = vec![vec![Rq::zero(); self.constants.KAPPA]; self.constants.R];

        if is_verbose() {
            println!("Generating inner Ajtai commitments...");
        }
        // Compute inner Ajtai commitments
        for kappa_iter in 0..self.constants.KAPPA {
            let a_mat_row = crs.fetch_next_n(self.constants.N);

            for i in 0..self.constants.R {
                let s_i = self.witness.column(i).to_vec();
                let t_i_elem = polynomial_vec_inner_product(&a_mat_row, &s_i);
                t_i_all[i][kappa_iter] = t_i_elem;
            }
        }
        
        if is_verbose() {
            println!("Computed Inner Ajtai Commitments!");
        }

        if is_verbose() {
            println!("Generating g_mat...");
        }
        // NOTE: This is a symmetric matrix
        let mut g_mat = Array2::from_elem((self.constants.R, self.constants.R), Rq::new(vec![]));
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                let col_i = self.witness.column(i).to_vec();
                let col_j = self.witness.column(j).to_vec();

                let res: Rq = polynomial_vec_inner_product(&col_i, &col_j);
                if g_mat[[i, j]] == Rq::new(vec![]) {
                    g_mat[[i, j]] = res
                }
            }
        }
        if is_verbose() {
            println!("Finished generating g_mat.");
        }

        if is_verbose() {
            println!("Generating polynomial decomposition of t_i");
        }
        let mut lhs = gen_empty_poly_vec(self.constants.KAPPA_1 as usize);
        for i in 0..self.constants.R {
            println!("i = {}", &i);
            let t_i_decomposed: Vec<Vec<Rq>> = decompose_polynomial_vec(&t_i_all[i], self.constants.B_1, self.constants.T_1);
            for k in 0..(self.constants.T_1 as usize) {
                println!("k = {}", &k);

                let mut prod : Vec<Rq> = vec![];
                
                for kappa_iter in 0..self.constants.KAPPA_1 {
                    let b_ik_row = crs.fetch_next_n(self.constants.KAPPA);
                    prod.push(polynomial_vec_inner_product(&b_ik_row, &t_i_decomposed[k]))
                }
                
                lhs = add_poly_vec(&prod, &lhs);
            }
        }

        if is_verbose() {
            println!("Generating polynomial decomposition of g_i");
        }
        let mut rhs = gen_empty_poly_vec(self.constants.KAPPA_2 as usize);
        for i in 0..self.constants.R {
            for j in i..self.constants.R {
                let g_ij: Vec<Rq> = decompose_polynomial(&g_mat[[i, j]], self.constants.B_2, self.constants.T_2);
                for k in 0..(self.constants.T_2 as usize) {
                    let c_ijk = crs.fetch_next_n(self.constants.KAPPA_2);
                    let g_ij_k: &Rq = &g_ij[k];
                    let prod = poly_by_poly_vec(g_ij_k, &c_ijk);
                    rhs = add_poly_vec(&prod, &rhs);
                }
            }
        }

        if is_verbose() {
            println!("Generating u_1...");
        }
        let u_1: Vec<Rq> = add_poly_vec(&lhs, &rhs);
        //println!("IMPORTANT PRINT!!! lhs: {:?}, rhs: {:?}", &lhs, &rhs);

        if is_verbose() {
            println!("Computing JL projection...");
        }
        // Next, compute JL projection

        let mut jl_tuple = self.jl_project();

        let mut projection_ints: Vec<i128> = jl_tuple.0;
        let mut pi_i_all: Vec<Array2<Zq>> = jl_tuple.1;

        if is_verbose() {
            println!("Computed projection..");
        }
        let mut rejections = 0;
        while !self.verifier.valid_projection(&projection_ints) {
            if is_verbose() {
                println!("Verifier rejected projection!");
            }
            rejections += 1;
            if rejections > 5 {
                panic!("failed JL...");
            }
            jl_tuple = self.jl_project();
            projection_ints = jl_tuple.0;
            pi_i_all = jl_tuple.1;
        }
        if is_verbose() {
            println!("JL Projection complete and accepted");
        }

        let projection = Zq::lift(&projection_ints);

        // First Aggregation step
        let upper_bound: usize = std::cmp::min(K, (128.0f64 / (*Q as f64).log2()).ceil() as usize);

        let mut psi: Vec<Vec<Zq>> = vec![];
        let mut omega: Vec<Vec<Zq>> = vec![];

        // TODO these comments aren't precise enough in dimensions
        // vector containing all phi'' for k in 1..(upper_bound), i in 0..(R-1)
        // TODO also perhaps switch this to an actual Vec<Array2<Rq>> later.
        let mut phi_prime_prime_vec: Vec<Vec<Vec<Rq>>> = vec![];

        // vector containing all b'' for k in 1..(upper_bound)
        let mut b_prime_prime_vec: Vec<Rq> = vec![];

        if is_verbose() {
            println!("First aggregation step...");
        }
        for k in 1..(upper_bound + 1) {
            if is_verbose() {
                println!("k={}, upper_bound: {}", k, upper_bound);
            }
            let psi_k = self.verifier.generate_psi();
            let omega_k = self.verifier.generate_omega();

            let mut phi_prime_prime_k: Vec<Vec<Rq>> = vec![];

            psi.push(psi_k.clone());
            omega.push(omega_k);

            // NOTE: for our first basic rudimentary implementation, consider that
            // a_prime_ij is just going to be a_ij... given that our F' is just going to be F since
            // they both have zero constant coefficient. Simple F' for now.

            // TODO obviously we'll iterate the below look across all L in the end,
            // AND a_prime_prime will have a (k) superscript... so this will be a full K length vec eventually
            let mut a_prime_prime: Array2<Rq> = Array2::zeros((self.constants.R, self.constants.R));

            for i in 0..self.constants.R {
                let phi_i_prime_prime: Vec<Rq>;
                for j in 0..self.constants.R {
                    let a_prime_ij: &Rq = &st.a_prime_k[0][[i, j]];
                    let a_prime_prime_ij: Rq = multiply_poly_ints(&a_prime_ij, &psi[k - 1]);
                    a_prime_prime[[i, j]] = a_prime_prime_ij;
                }

                // TODO we're taking the k=0 index for now.. disregarding the rest
                let phi_prime_i: Vec<Rq> = st.phi_k[0].column(i).to_vec();

                let mut lhs: Vec<Rq> = gen_empty_poly_vec(self.constants.N);
                for _l in 0..L {
                    lhs = add_poly_vec(
                        &multiply_poly_vec_ints(&phi_prime_i, &psi.last().unwrap()),
                        &lhs,
                    );
                }

                // side note: consider that pi_i^(j) is the jth row of pi_i for j = 1, ..., 256.
                let pi_i = &pi_i_all[i];
                let mut rhs: Vec<Rq> = vec![Rq::new(vec![]); self.constants.N as usize];
                for j in 0..256 {
                    let bolded_pi_poly_vec: Vec<Rq> = concat_coeff_reduction(&pi_i.row(j).to_vec());
                    let conj = sigma_inv_vec(&bolded_pi_poly_vec);
                    let omega_k_j = &omega[0][j];
                    let res = scale_poly_vec(&conj, omega_k_j);
                    rhs = add_poly_vec(&rhs, &res);
                }
                phi_i_prime_prime = add_poly_vec(&lhs, &rhs);
                phi_prime_prime_k.push(phi_i_prime_prime);
            }

            let mut lhs = Rq::new(vec![Zq::zero()]);
            for i in 0..self.constants.R {
                for j in 0..self.constants.R {
                    let prod = polynomial_vec_inner_product(
                        &self.witness.column(i).to_vec(),
                        &self.witness.column(j).to_vec(),
                    );
                    let res = &a_prime_prime[[i, j]] * prod;
                    lhs = lhs + res;
                }
            }

            let mut rhs = Rq::new(vec![Zq::zero()]);
            for i in 0..self.constants.R {
                let res = polynomial_vec_inner_product(
                    &phi_prime_prime_k[i],
                    &self.witness.column(i).to_vec(),
                );
                rhs = rhs + res;
            }
            let b_prime_prime_k = lhs.clone() + rhs.clone();

            self.verifier.verify_b_prime_prime(
                &b_prime_prime_k,
                &omega.last().unwrap(),
                &psi.last().unwrap(),
                &projection,
            );

            phi_prime_prime_vec.push(phi_prime_prime_k);
            b_prime_prime_vec.push(b_prime_prime_k);
        }

        if is_verbose() {
            println!("fetching alpha and beta...");
        }

        let alpha: Vec<Rq> = self.verifier.fetch_alpha();
        let beta: Vec<Rq> = self.verifier.fetch_beta();

        if is_verbose() {
            println!("Generating phi final..");
        }
        let mut phi_final: Vec<Vec<Rq>> = vec![];
        for i in 0..self.constants.R {
            let mut lhs: Vec<Rq> = gen_empty_poly_vec(self.constants.N);
            for k in 0..K {
                let res = poly_by_poly_vec(&alpha[k], &st.phi_k[k].column(i).to_vec());
                lhs = add_poly_vec(&lhs, &res);
            }
            let mut rhs: Vec<Rq> = gen_empty_poly_vec(self.constants.N);
            for k in 0..upper_bound {
                let res = poly_by_poly_vec(&beta[k], &phi_prime_prime_vec[k][i]);
                rhs = add_poly_vec(&rhs, &res);
            }
            phi_final.push(add_poly_vec(&lhs, &rhs));
        }

        if is_verbose() {
            println!("Generating h_mat..");
        }
        // NOTE: This is *also* a symmetric matrix
        let mut h_mat = Array2::from_elem((self.constants.R, self.constants.R), Rq::new(vec![]));
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                //println!("i : {}, j: {}", i, j);
                let s_i = self.witness.column(i).to_vec();
                let s_j = self.witness.column(j).to_vec();
                //println!("s_i: {:?}, s_j: {:?}", s_i, s_j);

                let phi_i = &phi_final[i];
                let phi_j = &phi_final[j];
                //println!("phi_i: {:?}, phi_j: {:?}", phi_i, phi_j);

                MOD_SUSPENSION.store(true, Ordering::SeqCst);

                let sum = polynomial_vec_inner_product(phi_i, &s_j)
                    + polynomial_vec_inner_product(phi_j, &s_i);

                //println!("sum: {}", sum);

                // NOTE here we're computing the
                // multiplicative inverse of 2, i.e., "1/2"
                let two = BigInt::from(2);
                let q_bigint = BigInt::from(*Q);
                let inv: i128 = two
                    .modpow(&(q_bigint.clone() - BigInt::from(2)), &q_bigint)
                    .to_i128()
                    .unwrap();
                let mut res = scale_polynomial(&sum, &Zq::from(inv));

                MOD_SUSPENSION.store(false, Ordering::SeqCst);
                res = res.recompute_mod();

                //println!("res: {}", res);

                if h_mat[[i, j]] == Rq::new(vec![]) {
                    h_mat[[i, j]] = res;
                }
            }
        }

        if is_verbose() {
            println!("Computing u_2...");
        }
        // Compute u_2
        let mut u_2: Vec<Rq> = vec![Rq::new(vec![]); self.constants.KAPPA_2 as usize];
        for i in 0..self.constants.R {
            for j in i..self.constants.R {
                for k in 0..(self.constants.T_1 as usize) {
                    let d_ijk_vec = crs.fetch_next_n(self.constants.KAPPA_2);
                    // TODO I don't think we can just add as such.
                    let dec_hij = decompose_polynomial(&h_mat[[i, j]], self.constants.B_1, self.constants.T_1);
                    let prod = poly_by_poly_vec(&dec_hij[k], &d_ijk_vec);
                    // NOTE prod.len() = KAPPA_2 (it should at least)
                    for idx in 0..prod.len() {
                        u_2[idx] = &u_2[idx] + &prod[idx];
                    }
                }
            }
        }

        let mut z: Vec<Rq> = vec![Rq::new(vec![]); self.constants.N as usize];
        let mut c_vec: Vec<Rq> = vec![];

        if is_verbose() {
            println!("challenge polynomial fetching..");
        }
        //println!("JUST FOR SANITY, WITNESS SHAPE: {:?}", &self.S.shape());
        for i in 0..self.constants.R {
            let c_i = self.verifier.fetch_challenge();
            c_vec.push(c_i);
            let prod: Vec<Rq> = poly_by_poly_vec(
                c_vec.last().as_ref().unwrap(),
                &self.witness.column(i).to_vec(),
            );
            //println!("WITNESS COL VEC SIZE: {}", &self.S.column(i).to_vec().len());
            //println!("PROD VEC SIZE: {}", &prod.len());
            for j in 0..self.constants.N {
                z[j] = &z[j] + &prod[j];
            }
        }
        //println!("z length.. {}" , &z.len());
        if is_verbose() {
            println!("Filling proof transcript.");
        }


        crs.reset_offset(); // we're about to verify the proof, so we need to be able to
                            // reconstruct the the full CRS deterministically from the base seed.

        // TODO fill the proof transcript with all the relevant data and return
        let proof_transcript: Transcript = Transcript {
            u_1,
            pi_i_all,
            projection,
            psi,
            omega,
            b_prime_prime: b_prime_prime_vec,
            alpha,
            beta,
            u_2,
            c: c_vec,
            z,
            t_i_all,
            g_mat,
            h_mat,
        };
        proof_transcript
    }

    pub fn jl_project(&mut self) -> (Vec<i128>, Vec<Array2<Zq>>) {
        // verifier sends random matrices in {-1,0,1}...
        let mut projection: Vec<i128> = vec![0; 256];
        let mut pi_i_all: Vec<Array2<Zq>> = vec![];
        for i in 0..self.constants.R {
            let pi_i: Array2<i128> = self.verifier.sample_jl_projection();
            let s_i_coeffs: Array2<i128> = vec_to_column_array(&Zq::lift_inv(
                &witness_coeff_concat(&self.witness.column(i).to_vec()),
            ));

            // NOTE: for reference, this is a 256x(ND) multiplied by an (ND)x1, giving a 256x1
            // which we turn into a vec

            let product = matmul(&pi_i, &s_i_coeffs).column(0).to_vec();
            projection = add_vecs(&projection, &product);

            let mut pi_i_lifted = Array2::from_elem((256, self.constants.N * (D as usize)), Zq::zero());
            // we have to lift Pi_i to make it suitable for computations in Rq / Zq..
            for j in 0..256 {
                pi_i_lifted
                    .row_mut(j)
                    .assign(&ArrayView::from(&Zq::lift(&pi_i.row(j).to_vec())));
            }

            pi_i_all.push(pi_i_lifted);
        }
        (projection, pi_i_all)
    }
}

// generates and returns the MSIS witness S (here, just denoted 'witness')
pub fn generate_witness(constants: &RuntimeConstants) -> Array2<Rq> {
    let mut witness = Array2::from_elem(Ix2(constants.N, constants.R), Rq::new(vec![]));
    let mut rng = rand::thread_rng();

    // Random Generation of witness matrix

    // total summed norm:
    let mut norm_sum: i128 = 0;

    for i in 0..constants.N {
        for j in 0..constants.R {
            let witness_ij = generate_polynomial(*Q, D);
            // TODO fix clone
            witness[[i, j]] = witness_ij.clone();
            // add norm to norm sum
            norm_sum += poly_norm(&witness_ij) as i128;
        }
    }

    // if too big, scale down:
    while norm_sum > (i128::pow(constants.BETA_BOUND, 2)) {
        if is_verbose() {
            println!(
                "norm sum too big! {} {}",
                norm_sum,
                (constants.BETA_BOUND) * (constants.BETA_BOUND)
            );
        }

        let i = rng.gen_range(0..constants.N);
        let j = rng.gen_range(0..constants.R);

        let old_poly_norm = poly_norm(&witness[[i, j]]) as i128;

        //println!("i {}, j {},.... {}", i, j, poly_norm(&S[[i,j]]) as i128);
        //println!("POLY NORM: {}", poly_norm(&S[[i,j]]) as i128);
        //println!("norm sum: {}", norm_sum);
        //println!("Old Sij: {}", &S[[i,j]]);

        let witness_ij = reduce_polynomial(&witness[[i, j]]);

        //println!("New Sij: {}", &s_ij);
        let new_poly_norm = poly_norm(&witness_ij) as i128;
        //println!("norm diff... {}", old_poly_norm - new_poly_norm);
        //let mut s_ij = Rq::zero();

        witness[[i, j]] = witness_ij.clone();

        norm_sum -= old_poly_norm - new_poly_norm;
    }
    if is_verbose() {
        println!(
            "norm of witness! squared version: {}, sqrt version: {}",
            norm_sum,
            (norm_sum as f64).sqrt() as i128
        );
    }
    witness
}
