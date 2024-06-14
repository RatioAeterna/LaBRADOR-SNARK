use ndarray::Array2;
use rand::distributions::{Distribution, WeightedIndex};
use rand::Rng;

use crate::algebraic::*;
use crate::constants::*;
use crate::structs::*;
use crate::util::*;
use num_traits::Zero;
use rayon::prelude::*;

pub struct Verifier<'a> {
    b_prime: Option<Vec<Zq>>,
    constants: &'a RuntimeConstants,
}

impl<'a> Verifier<'a> {
    pub fn new(b_prime_k: Vec<Zq>, constants: &'a RuntimeConstants) -> Self {
        Verifier {
            b_prime: Some(b_prime_k),
            constants: constants,
        }
    }

    pub fn verify(&self, st: &State, proof: &Transcript, crs: &mut CRS) -> bool {
        // LINE 1, LINE 2
        // Enumeration of the state and transcript, respectively, so we don't include here.
        let upper_bound: usize = std::cmp::min(K, (128.0f64 / (*Q as f64).log2()).ceil() as usize);
        if is_verbose() {
            println!("Lines 1 and 2... enumeration of the state and the proof transcript");
        }

        if is_verbose() {
            println!("starting line 3");
        }
        // LINE 3
        // Computing "a_prime_prime" matrices for all k up to upper_bound, storing those in a vector
        let mut a_prime_prime: Vec<Array2<Rq>> = vec![];
        for k in 0..upper_bound {
            let mut a_prime_prime_mat = Array2::from_elem((self.constants.R, self.constants.R), Rq::new(vec![]));
            for i in 0..self.constants.R {
                for j in 0..self.constants.R {
                    let mut sum: Rq = Rq::new(vec![]);
                    for l in 0..L {
                        let scaled: Rq =
                            scale_polynomial(&st.a_prime_k[l][[i, j]], &proof.psi[k][l]);
                        sum = sum + scaled;
                    }
                    a_prime_prime_mat[[i, j]] = sum;
                }
            }
            a_prime_prime.push(a_prime_prime_mat);
        }

        if is_verbose() {
            println!("starting line 4");
        }
        // LINE 4
        // Computing "phi_i_prime_prime" vecs for each k (and all i's)
        let mut phi_prime_prime_k: Vec<Vec<Vec<Rq>>> = vec![];
        for k in 0..upper_bound {
            //println!("k={}", k);
            // contains all phi_prime_prime_i for this particular k
            let mut phi_prime_prime: Vec<Vec<Rq>> = vec![];
            for i in 0..self.constants.R {
                // ith polynomial vec for this particular k
                let mut sum: Vec<Rq> = vec![];
                for l in 0..L {
                    // TODO yes, we'll make this faster eventually
                    //println!("k: {}, l: {}, phi_prime_k len: {}, psi len: {}, psi_0 len: {}", k, l, &st.phi_prime_k.len(), proof.psi.len(), proof.psi[0].len());
                    let prod =
                        scale_poly_vec(&st.phi_prime_k[l].column(i).to_vec(), &proof.psi[k][l]);
                    if sum.len() == 0 {
                        sum = vec![Rq::zero(); prod.len()];
                    }

                    sum = add_poly_vec(&sum, &prod);
                }
                for j in 0..256 {
                    let bolded_pi_poly_vec: Vec<Rq> =
                        concat_coeff_reduction(&proof.pi_i_all[i].row(j).to_vec());
                    let conj = sigma_inv_vec(&bolded_pi_poly_vec);
                    let prod = scale_poly_vec(&conj, &proof.omega[k][j]);
                    sum = add_poly_vec(&sum, &prod);
                }
                phi_prime_prime.push(sum);
            }
            phi_prime_prime_k.push(phi_prime_prime);
        }

        if is_verbose() {
            println!("starting line 5");
        }
        // LINE 5
        // Forming a single canonical matrix of "a_ij" polynomials
        let mut a_constraints = Array2::from_elem((self.constants.R, self.constants.R), Rq::new(vec![]));
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                // generate single a_ij
                let mut a_ij: Rq = Rq::new(vec![]);
                for k in 0..K {
                    let prod: Rq = &proof.alpha[k] * &st.a_k[k][[i, j]];
                    a_ij = a_ij + prod;
                }
                for k in 0..upper_bound {
                    let prod: Rq = &proof.beta[k] * &a_prime_prime[k][[i, j]];
                    a_ij = a_ij + prod;
                }
                a_constraints[[i, j]] = a_ij;
            }
        }

        if is_verbose() {
            println!("starting line 6");
        }
        // LINE 6
        // Forming a single canonical set of vectors "phi_i" for i in {0, ..., R-1}
        let mut phi: Vec<Vec<Rq>> = vec![];
        for i in 0..self.constants.R {
            let mut phi_i: Vec<Rq> = vec![];
            for k in 0..K {
                let prod = poly_by_poly_vec(&proof.alpha[k], &st.phi_k[k].column(i).to_vec());
                if phi_i.len() == 0 {
                    phi_i = gen_empty_poly_vec(prod.len());
                }
                phi_i = add_poly_vec(&phi_i, &prod);
            }
            for k in 0..upper_bound {
                let prod = poly_by_poly_vec(&proof.beta[k], &phi_prime_prime_k[k][i]);
                phi_i = add_poly_vec(&phi_i, &prod);
            }
            phi.push(phi_i);
        }

        if is_verbose() {
            println!("starting line 7");
        }
        // LINE 7
        // Forming a single canonical "b" polynomial
        let mut b: Rq = Rq::new(vec![]);
        for k in 0..K {
            let prod: Rq = &proof.alpha[k] * &st.b_k[k];
            b = b + prod;
        }
        for k in 0..upper_bound {
            let prod: Rq = &proof.beta[k] * &proof.b_prime_prime[k];
            b = b + prod;
        }

        if is_verbose() {
            println!("starting line 8");
        }
        // CHECK 8
        // check that g_{ij} == g_{ji} i.e., matrix g_mat is symmetric
        // TODO is it faster to only check some values? I really don't think this makes a
        // difference.
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                if proof.g_mat[[i, j]] != proof.g_mat[[j, i]] {
                    return false;
                }
            }
        }

        if is_verbose() {
            println!("starting line 9");
        }
        // CHECK 9
        // check that h_{ij} == h{ji} i.e., matrix h_mat is symmetric
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                if proof.h_mat[[i, j]] != proof.h_mat[[j, i]] {
                    //println!("At i,j: {}", proof.h_mat[[i,j]]);
                    //println!("At j,i: {}", proof.h_mat[[j,i]]);
                    return false;
                }
            }
        }

        if is_verbose() {
            println!("starting line 10");
        }
        // LINE 10
        // Decompose vec z into z = z^(0) + z^(1)b
        let z_decompositions: Vec<Vec<Rq>> = decompose_polynomial_vec(&proof.z, self.constants.B, 2);

        if is_verbose() {
            println!("starting line 11");
        }
        // LINE 11
        // Decompose vec t_i the same way
        let mut t_decompositions: Vec<Vec<Vec<Rq>>> = vec![];
        for i in 0..self.constants.R {
            let t_i_decomposed: Vec<Vec<Rq>> =
                decompose_polynomial_vec(&proof.t_i_all[i], self.constants.B_1, self.constants.T_1);
            t_decompositions.push(t_i_decomposed);
        }

        if is_verbose() {
            println!("starting line 12");
        }
        // LINE 12
        // Decompose matrix elements g_ij
        let mut g_mat_decompositions = Array2::from_elem((self.constants.R, self.constants.R), vec![]);
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                //println!("g_mat term at [i={},j={}]... {}", i,j,&proof.g_mat[[i,j]]);
                let dec_gij: Vec<Rq> = decompose_polynomial(&proof.g_mat[[i, j]], self.constants.B_2, self.constants.T_2);
                g_mat_decompositions[[i, j]] = dec_gij;
                //println!("decomposition.... {:?}", g_mat_decompositions[[i,j]]);
            }
        }

        if is_verbose() {
            println!("starting line 13");
        }
        // LINE 13
        // Decompose matrix elements h_ij
        let mut h_mat_decompositions = Array2::from_elem((self.constants.R, self.constants.R), vec![]);
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                let dec_hij: Vec<Rq> = decompose_polynomial(&proof.h_mat[[i, j]], self.constants.B_1, self.constants.T_1);
                h_mat_decompositions[[i, j]] = dec_hij;
            }
        }

        if is_verbose() {
            println!("starting line 14");
        }
        // LINE 14
        // TODO Yes, we can flatten a number of these loops. Just want to get the protocol down
        // now, will save on computation later.
        let mut sum: f64 = 0.0;
        for i in 0..2 {
            sum += vec_poly_norm_squared(&z_decompositions[i]);
        }
        for i in 0..self.constants.R {
            for k in 0..(self.constants.T_1 as usize) {
                sum += vec_poly_norm_squared(&t_decompositions[i][k]);
            }
        }
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                //println!("g_mat term at [i={},j={}]... {}", i,j,&proof.g_mat[[i,j]]);
                for k in 0..(self.constants.T_2 as usize) {
                    sum += poly_norm(&g_mat_decompositions[[i, j]][k]);
                    //println!("g_mat decomposition at k={}, is: {}", k,&g_mat_decompositions[[i,j]][k]);
                }
            }
        }
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                for k in 0..(self.constants.T_1 as usize) {
                    sum += poly_norm(&h_mat_decompositions[[i, j]][k]);
                }
            }
        }
        //println!("sum: {}, should be less than... {}", sum, (*BETA_PRIME).powi(2));

        if sum > (self.constants.BETA_PRIME).powi(2) {
            return false;
        }

        if is_verbose() {
            println!("starting line 15");
        }

        // CHECK 15
        let mut lhs : Vec<Rq> = vec![];
        for kappa_iter in 0..self.constants.KAPPA {
            let a_mat_row = crs.fetch_A_row(kappa_iter);
            //println!("Here's an 'a' row! {:?}", &a_mat_row);
            lhs.push(polynomial_vec_inner_product(&a_mat_row, &proof.z));
        }
        let mut rhs: Vec<Rq> = vec![];
        for i in 0..self.constants.R {
            let prod = poly_by_poly_vec(&proof.c[i], &proof.t_i_all[i]);
            if rhs.len() == 0 {
                rhs = gen_empty_poly_vec(prod.len());
            }
            rhs = add_poly_vec(&rhs, &prod);
            //println!("c_i: {}", &proof.c[i]);
            //println!("t_i: {:?}", &proof.t_i_all[i]);
        }

        if lhs != rhs {
            //println!("LHS: {:?}", lhs);
            //println!("RHS: {:?}", rhs);
            //println!("z: {:?}", &proof.z);
            return false;
        }

        if is_verbose() {
            println!("starting line 16");
        }

        // CHECK 16
        let lhs: Rq = polynomial_vec_inner_product(&proof.z, &proof.z);
        let mut rhs: Rq = Rq::new(vec![]);

        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                rhs = rhs + (&proof.g_mat[[i, j]] * &proof.c[i] * &proof.c[j]);
            }
        }

        if lhs != rhs {
            return false;
        }

        if is_verbose() {
            println!("starting line 17");
        }
        // CHECK 17
        let mut lhs: Rq = Rq::new(vec![]);
        //println!("{} {}", &proof.c.len(), &phi.len());
        for i in 0..self.constants.R {
            lhs = lhs + polynomial_vec_inner_product(&phi[i], &proof.z) * &proof.c[i];
        }
        let mut rhs: Rq = Rq::new(vec![]);
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                rhs = rhs + (&proof.h_mat[[i, j]] * &proof.c[i] * &proof.c[j]);
            }
        }

        if lhs != rhs {
            return false;
        }

        if is_verbose() {
            println!("starting line 18");
        }
        // CHECK 18
        let mut s1: Rq = Rq::new(vec![]);
        let mut s2: Rq = Rq::new(vec![]);
        for i in 0..self.constants.R {
            for j in 0..self.constants.R {
                let prod: Rq = &a_constraints[[i, j]] * &proof.g_mat[[i, j]];
                s1 = s1 + prod;
            }
            s2 = s2 + &proof.h_mat[[i, i]];
        }
        // check to make sure this is the zero polynomial
        if (s1 + s2 - b) != Rq::new(vec![]) {
            return false;
        }

        if is_verbose() {
            println!("starting line 19");
        }
        // CHECK 19
        //let mut u_1_candidate : Vec<Rq> = vec![Rq::new(vec![]); KAPPA_1 as usize];
        /*
        let mut lhs = gen_empty_poly_vec(self.constants.KAPPA_1 as usize);
        for i in 0..self.constants.R {
            for k in 0..(self.constants.T_1 as usize) {
                let mut prod : Vec<Rq> = vec![];
                for kappa_iter in 0..self.constants.KAPPA_1 {
                    let b_ik_row = crs.fetch_B_ik_row(i, k, kappa_iter);
                    prod.push(polynomial_vec_inner_product(&b_ik_row, &t_decompositions[i][k]))
                }
                lhs = add_poly_vec(&prod, &lhs);
            }
        }
        */
        let intermediate_results: Vec<Vec<Vec<Rq>>> = (0..self.constants.R)
            .into_par_iter()
            .map(|i| {
                (0..self.constants.T_1 as usize)
                    .into_par_iter()
                    .map(|k| {
                        let prod: Vec<Rq> = (0..self.constants.KAPPA_1)
                            .map(|kappa_iter| {
                                let b_ik_row = crs.fetch_B_ik_row(i, k, kappa_iter);
                                polynomial_vec_inner_product(&b_ik_row, &t_decompositions[i][k])
                            })
                            .collect();
                        prod
                    })
                    .collect()
            })
            .collect();

        // Sequentially combine the results
        let mut lhs = gen_empty_poly_vec(self.constants.KAPPA_1 as usize);
        for results in intermediate_results {
            for prod in results {
                lhs = add_poly_vec(&prod, &lhs);
            }
        }

        let mut rhs = gen_empty_poly_vec(self.constants.KAPPA_2 as usize);
        for i in 0..self.constants.R {
            for j in i..self.constants.R {
                for k in 0..(self.constants.T_2 as usize) {
                    let c_ijk = crs.fetch_C_ijk(i, j, k);
                    let poly = &g_mat_decompositions[[i, j]][k];
                    let prod = poly_by_poly_vec(poly, &c_ijk);
                    rhs = add_poly_vec(&prod, &rhs);
                }
            }
        }
        let u_1_candidate: Vec<Rq> = add_poly_vec(&lhs, &rhs);

        //println!("PROOF u1: {:?}, COMPUTED u1... {:?}", &proof.u_1, &u_1_candidate);
        //println!("lhs: {:?}, rhs: {:?}", &lhs, &rhs);
        if &proof.u_1 != &u_1_candidate {
            return false;
        }

        if is_verbose() {
            println!("starting line 20");
        }
        // CHECK 20
        let mut u_2_candidate: Vec<Rq> = vec![Rq::new(vec![]); self.constants.KAPPA_2];
        for i in 0..self.constants.R {
            for j in i..self.constants.R {
                for k in 0..(self.constants.T_1 as usize) {
                    let d_ijk_vec = crs.fetch_D_ijk(i,j,k);
                    let prod = poly_by_poly_vec(&h_mat_decompositions[[i, j]][k], &d_ijk_vec);
                    // NOTE prod.len() = KAPPA_2 (it should at least)
                    u_2_candidate = add_poly_vec(&u_2_candidate, &prod);
                }
            }
        }
        // now, check for equality with actual u_2
        if &proof.u_2 != &u_2_candidate {
            return false;
        }
        //crs.reset_offset(); // reset the CRS offset seed in case it is re-used later
        true
    }

    // TODO do we want to STORE alpha, beta in the Verifier struct?
    pub fn fetch_alpha(&self) -> Vec<Rq> {
        let mut alpha = vec![];
        for _i in 0..K {
            alpha.push(generate_polynomial(*Q, D));
        }
        alpha
    }

    pub fn fetch_beta(&self) -> Vec<Rq> {
        let mut beta = vec![];
        let upper_bound: usize = std::cmp::min(K, (128.0f64 / (*Q as f64).log2()).ceil() as usize);
        for _i in 0..upper_bound {
            beta.push(generate_polynomial(*Q, D));
        }
        beta
    }

    // fetch a challenge polynomial from the challenge space \mathcal{C} satisfying a number of
    // criteria
    pub fn fetch_challenge(&self) -> Rq {
        // particular challenge coefficient distribution described on page 6.
        let mut coeff_dist: Vec<Zq> = vec![];
        if D == 64 {
            for _i in 0..23 {
                coeff_dist.push(Zq::from(0));
            }
            for _i in 0..31 {
                coeff_dist.push(Zq::from(1));
            }
            for _i in 0..10 {
                coeff_dist.push(Zq::from(2));
            }
        } else {
            coeff_dist.push(Zq::from(1));
            coeff_dist.push(Zq::from(0));
            //coeff_dist.push(Zq::from(0));
        }

        let mut candidate = generate_polynomial_picky(D as usize, coeff_dist.clone());
        // TODO... I don't think this norm should be squared, as it is.. which would perhaps give you 71^2.. definitely fix this if
        // needed.
        //assert!(poly_norm(&candidate) == TAU, "Incorrect l2 norm of challenge polynomial");

        while operator_norm(&candidate) > T {
            //assert!(poly_norm(&candidate) == TAU, "Incorrect l2 norm of challenge polynomial");
            candidate = generate_polynomial_picky(D as usize, coeff_dist.clone());
        }
        candidate
    }

    pub fn generate_psi(&self) -> Vec<Zq> {
        let mut rng = rand::thread_rng();
        let mut psi_k: Vec<Zq> = Vec::new();
        for _i in 0..L {
            psi_k.push(Zq::from(rng.gen_range(0..*Q)));
            //println!("PSI K: {:?}", &psi_k);
            // TODO for debug purposes
            //psi_k.push(Zq::from(0));
        }
        psi_k
    }

    pub fn generate_omega(&self) -> Vec<Zq> {
        let mut rng = rand::thread_rng();
        let mut omega_k: Vec<Zq> = Vec::new();
        for _i in 0..256 {
            omega_k.push(Zq::from(rng.gen_range(0..*Q)));
            // TODO for debug purposes
            //omega_k.push(Zq::from(0));
        }
        //self.omega_k = Some(&omega_k);
        omega_k
    }

    pub fn fetch_alleged_b_prime_prime_cc(
        &self,
        omega_k: &Vec<Zq>,
        psi_k: &Vec<Zq>,
        projection: &Vec<Zq>,
    ) -> Zq {
        //let prod = Zq::from(vec_inner_product(&Zq::lift_inv(omega_k), projection));
        let prod = Zq::from(vec_inner_product_zq(omega_k, projection));
        let mut sum = Zq::zero();
        for i in 0..L {
            sum += &psi_k[i] * &(self.b_prime.as_ref().unwrap()[i]);
        }
        let alleged_b_prime_prime = prod + sum;

        alleged_b_prime_prime
    }

    pub fn verify_b_prime_prime(
        &self,
        b_prime_prime_k: &Rq,
        omega_k: &Vec<Zq>,
        psi_k: &Vec<Zq>,
        projection: &Vec<Zq>,
    ) -> () {
        // TODO again column vs row not sure.
        // Also self vs no self keyword not sure.
        //let prod = Zq::from(vec_inner_product(&Zq::lift_inv(omega_k), projection));
        let prod = Zq::from(vec_inner_product_zq(omega_k, projection));
        let mut sum = Zq::zero();
        for i in 0..L {
            sum += &psi_k[i] * &(self.b_prime.as_ref().unwrap()[i]);
        }
        let check_val = prod + sum;

        // check that the constant term is equal to the above stuff.
        assert!(b_prime_prime_k.eval(Zq::zero()) == check_val, "verify_b_prime_prime check failed: b_prime_prime_k constant coefficient is: {}, and check_val is: {}", b_prime_prime_k.eval(Zq::zero()), check_val);
    }

    pub fn sample_jl_projection(&self) -> Array2<i128> {
        let mut rng = rand::thread_rng();
        let choices = [-1, 0, 1];
        let weights = [0.25, 0.5, 0.25];
        let dist = WeightedIndex::new(&weights).unwrap();

        let mut pi_i: Array2<i128> = Array2::zeros((256, self.constants.N * (D as usize)));

        for ((_i, _j), value) in pi_i.indexed_iter_mut() {
            *value = choices[dist.sample(&mut rng)] as i128;
        }

        pi_i
    }

    pub fn valid_projection(&self, projection: &Vec<i128>) -> bool {
        let val: f64 = 128.;
        let norm = l2_norm(projection);
        if is_verbose() {
            println!(
                "TOTAL NORM OF JL PROJECTION: {}, {}",
                norm,
                val.sqrt() * (self.constants.BETA_BOUND as f64)
            );
        }
        return norm <= val.sqrt() * (self.constants.BETA_BOUND as f64);
    }
}
