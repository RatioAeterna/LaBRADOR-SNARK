use ndarray::{Array2, Ix2, concatenate};
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::{Distribution, WeightedIndex, Uniform};

use num_traits::Zero;
use crate::util::*;
use crate::algebraic::*;
use crate::constants::*;
use crate::structs::*;

pub struct Verifier {

    //projection : Option<&'a Array2<R_q>>,
    //psi_k : Option<Vec<Z_q>>,
    //omega_k : Option<&'a Vec<Z_q>>,
    b_prime : Option<Vec<Z_q>>,

    // JL projection matrix
    // elements of i128, not Z_q, since it's a bunch of random uniform {-1,0,1}
    Pi : Option<Vec<Array2<i128>>>,
}

impl Verifier {

    pub fn new(b_prime_k: Vec<Z_q>) -> Self {
        Verifier {
            //projection: None,
            //psi_k: None,
            //omega_k: None,
            b_prime: Some(b_prime_k),
            Pi: Some(vec![]),
        }
    }


    pub fn verify(&self, st: &State, proof : Transcript, crs : &CRS) -> bool {
        // LINE 1, LINE 2
        // Enumeration of the state and transcript, respectively, so we don't include here.
        let upper_bound : usize = (128.0f64 / (Q as f64).log10()).ceil() as usize;

        println!("starting line 3");
        // LINE 3
        // Computing "a_prime_prime" matrices for all k up to upper_bound, storing those in a vector

        let mut a_prime_prime : Vec<Array2<R_q>> = vec![];
        for k in 0..upper_bound {
            let mut a_prime_prime_mat = Array2::from_elem((R,R), R_q::new(vec![])); 
            for i in 0..R {
                for j in 0..R {
                    let mut sum : R_q = R_q::new(vec![]);
                    for l in 0..L {
                        let scaled : R_q = scale_polynomial(&st.a_prime_k[l][[i,j]], f32::from(proof.psi[k][l]));
                        sum = sum + scaled;
                    }
                    a_prime_prime_mat[[i,j]] = sum;
                }
            }
            a_prime_prime.push(a_prime_prime_mat);
        }

        println!("starting line 4");
        // LINE 4
        // Computing "phi_i_prime_prime" vecs for each k (and all i's)
        let mut phi_prime_prime_k : Vec<Vec<Vec<R_q>>> = vec![];
        for k in 0..upper_bound {
            // contains all phi_prime_prime_i for this particular k
            let mut phi_prime_prime : Vec<Vec<R_q>> = vec![];
            for i in 0..R {
                // ith polynomial vec for this particular k
                let mut sum : Vec<R_q> = vec![];
                for l in 0..L {
                    // TODO yes, we'll make this faster eventually
                    println!("k: {}, phi_prime_k len: {}, psi len: {}, psi_0 len: {}", k, &st.phi_prime_k.len(), proof.psi.len(), proof.psi[0].len());
                    let prod = scale_poly_vec(&st.phi_prime_k[k].column(i).to_vec(), f32::from(proof.psi[k][l]));
                    if (sum.len() == 0) {
                        sum = vec![R_q::zero(); prod.len()];
                    }

                    sum = add_poly_vec(&sum, &prod);
                }
                for j in 0..256 {
                    let bolded_pi_poly_vec : Vec<R_q> = concat_coeff_reduction(&Z_q::lift(&self.get_Pi_i(i).row(j).to_vec()));
                    let conj = sigma_inv_vec(&bolded_pi_poly_vec);
                    let prod = scale_poly_vec(&conj, f32::from(proof.omega[k][j]));
                    sum = add_poly_vec(&sum, &prod);
                }
                phi_prime_prime.push(sum);
            }
            phi_prime_prime_k.push(phi_prime_prime);
        }

        println!("starting line 5");
        // LINE 5
        // Forming a single canonical matrix of "a_ij" polynomials
        let mut Aij = Array2::from_elem((R,R), R_q::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                // generate single a_ij
                let mut a_ij : R_q = R_q::new(vec![]);
                for k in 0..K {
                    let prod : R_q = &proof.alpha[k] * &st.a_k[k][[i,j]];
                    a_ij = a_ij + prod;
                }
                for k in 0..upper_bound {
                    let prod : R_q = &proof.beta[k] * &a_prime_prime[k][[i,j]];
                    a_ij = a_ij + prod;
                }
                Aij[[i,j]] = a_ij;
            }
        }

        println!("starting line 6");
        // LINE 6 
        // Forming a single canonical set of vectors "phi_i" for i in {0, ..., R-1}
        let mut phi : Vec<Vec<R_q>> = vec![];
        for i in 0..R {
            let mut phi_i : Vec<R_q> = vec![]; 
            for k in 0..K {
                let prod = poly_by_poly_vec(&proof.alpha[k], &st.phi_k[k].column(i).to_vec());
                phi_i = add_poly_vec(&phi_i, &prod);
            }
            for k in 1..upper_bound {
                let prod = poly_by_poly_vec(&proof.beta[k], &phi_prime_prime_k[k][i]);
                phi_i = add_poly_vec(&phi_i, &prod);
            }
        }

        println!("starting line 7");
        // LINE 7
        // Forming a single canonical "b" polynomial
        let mut b : R_q = R_q::new(vec![]);
        for k in 0..K {
            let prod : R_q = &proof.alpha[k] * &st.b_k[k];
            b = b + prod;
        }
        for k in 0..upper_bound {
            let prod : R_q = &proof.beta[k] * &proof.b_prime_prime[k];
            b = b + prod;
        }
        
        println!("starting line 8");
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

        println!("starting line 9");
        // CHECK 9
        // check that h_{ij} == h{ji} i.e., matrix Hij is symmetric
        for i in 0..R {
            for j in 0..R {
                if (proof.Hij[[i,j]] != proof.Hij[[j,i]]) {
                    return false;
                }
            }
        }


        // LINE 10
        // Decompose vec z into z = z^(0) + z^(1)b
        let z_decompositions : Vec<Vec<R_q>> = decompose_polynomial_vec(&proof.z, *B, 2);

        // LINE 11
        // Decompose vec t_i the same way
        let mut t_decompositions: Vec<Vec<Vec<R_q>>> = vec![];
        for i in 0..R {
            let t_i_decomposed : Vec<Vec<R_q>> = decompose_polynomial_vec(&proof.t_i_all[i], *B_1, *T_1);
            t_decompositions.push(t_i_decomposed);
        } 
        
        // LINE 12
        // Decompose matrix elements g_ij
        let mut Gij_decompositions = Array2::from_elem((R,R), vec![]); 
        for i in 0..R {
            for j in 0..R {
                let dec_gij: Vec<R_q> = decompose_polynomial(&proof.Gij[[i,j]], *B_2, *T_2);
                Gij_decompositions[[i,j]] = dec_gij;
            }
        }

        // LINE 13
        // Decompose matrix elements h_ij
        let mut Hij_decompositions = Array2::from_elem((R,R), vec![]); 
        for i in 0..R {
            for j in 0..R {
                let dec_hij: Vec<R_q> = decompose_polynomial(&proof.Hij[[i,j]], *B_1, *T_1);
                Hij_decompositions[[i,j]] = dec_hij;
            }
        }


        // LINE 14
        // TODO Yes, we can flatten a number of these loops. Just want to get the protocol down
        // now, will save on computation later.
        let mut sum : f64 = 0.0; 
        for i in 0..2 {
            sum += vec_poly_norm_squared(&z_decompositions[i]); 
        }
        for i in 0..R {
            for k in 0..(*T_1 as usize) {
                sum += vec_poly_norm_squared(&t_decompositions[i][k]); 
            }
        }
        for i in 0..R {
            for j in 0..R {
                for k in 0..(*T_2 as usize) {
                    sum += poly_norm(&Gij_decompositions[[i,j]][k]); 
                }
            }
        }
        for i in 0..R {
            for j in 0..R {
                for k in 0..(*T_1 as usize) {
                    sum += poly_norm(&Hij_decompositions[[i,j]][k]); 
                }
            }
        }
        if (sum > (*BETA_PRIME).powi(2)) { return false; }


        // CHECK 15
        let mut lhs = polynomial_matrix_product(&crs.A_mat, &vec_to_column_array(&proof.z)).column(0).to_vec();
        let mut rhs : Vec<R_q> = vec![];
        for i in 0..R {
            let prod = poly_by_poly_vec(&proof.c[i], &proof.t_i_all[i]);
            let rhs = add_poly_vec(&rhs, &prod);
        }
        if (lhs != rhs) { return false;}

        // CHECK 16
        let lhs : R_q = polynomial_vec_inner_product(&proof.z, &proof.z);
        let mut rhs : R_q = R_q::new(vec![]);

        for i in 0..R {
            for j in 0..R {
                rhs = rhs + (&proof.Gij[[i,j]] * &proof.c[i] * &proof.c[j]);
            }
        }
        if (lhs != rhs) {
            return false;
        }

        // CHECK 17
        let mut lhs : R_q = R_q::new(vec![]);
        for i in 0..R {
            lhs = lhs + polynomial_vec_inner_product(&phi[i], &proof.z) * &proof.c[i];
        }
        let mut rhs : R_q = R_q::new(vec![]);
        for i in 0..R {
            for j in 0..R {
                rhs = rhs + (&proof.Hij[[i,j]] * &proof.c[i] * &proof.c[j]);
            }
        }
        if (lhs != rhs) { return false;}
    
        // CHECK 18
        let mut s1 : R_q = R_q::new(vec![]);
        let mut s2 : R_q = R_q::new(vec![]);
        for i in 0..R {
            for j in 0..R {
                let prod : R_q = &Aij[[i,j]] * &proof.Gij[[i,j]];
                s1 = s1 + prod;
            }
            s2 = s2 + &proof.Hij[[i,i]];
        }
        // check to make sure this is the zero polynomial
        if ((s1 + s2 - b) != R_q::new(vec![])) { return false; }

        // CHECK 19
        let mut u_1_candidate : Vec<R_q> = vec![R_q::new(vec![]); KAPPA_1 as usize];
        for i in 0..R {
            for k in 0..(*T_1 as usize) {
                let B_ik = crs.B_mat.get(&(i,k)).unwrap();
                let col_mat = vec_to_column_array(&proof.t_i_all[i]);
                let prod = polynomial_matrix_product(&B_ik, &col_mat).column(0).to_vec();
                u_1_candidate = add_poly_vec(&u_1_candidate, &prod);
            }
        }

        for i in 0..R {
            for j in i..R {
                for k in 0..(*T_2 as usize) {
                    let C_ijk = crs.C_mat.get(&(i,j,k)).unwrap().column(0).to_vec();
                    let poly = &Gij_decompositions[[i,j]][k];
                    let prod = poly_by_poly_vec(poly, &C_ijk);
                    u_1_candidate = add_poly_vec(&u_1_candidate, &prod);
                }
            }
        }
        if (&proof.u_1 != &u_1_candidate) { return false; }

        // CHECK 20
        let mut u_2_candidate : Vec<R_q> = vec![R_q::new(vec![]); KAPPA_2 as usize];
        for i in 0..R {
            for j in i..R {
                for k in 0..(*T_1 as usize) {
                    let D_ijk_vec = crs.D_mat.get(&(i,j,k)).unwrap().column(0).to_vec();
                    let prod = poly_by_poly_vec(&Hij_decompositions[[i,j]][k], &D_ijk_vec);
                    // NOTE prod.len() = KAPPA_2 (it should at least)
                    u_2_candidate = add_poly_vec(&u_2_candidate, &prod);
                }
            }
        }
        // now, check for equality with actual u_2
        if (&proof.u_2 != &u_2_candidate) { return false; }
        true
    }

    
    pub fn get_Pi_i(&self, i : usize) -> &Array2<i128> {
        &self.Pi.as_ref().unwrap()[i]
    }


    // TODO do we want to STORE alpha, beta in the Verifier struct?
    pub fn fetch_alpha(&self) -> Vec<R_q> {
        let mut alpha = vec![]; 
        for i in 0..K {
            alpha.push(generate_polynomial(Q,D));
        }
        alpha
    }

    pub fn fetch_beta(&self) -> Vec<R_q> {
        let mut beta = vec![]; 
        let upper_bound : usize = (128.0f64 / (Q as f64).log10()).ceil() as usize;
        for i in 0..upper_bound {
            beta.push(generate_polynomial(Q,D));
        }
        beta
    }


    // fetch a challenge polynomial from the challenge space \mathcal{C} satisfying a number of
    // criteria
    pub fn fetch_challenge(&self) -> R_q {
        // particular challenge coefficient distribution described on page 6.
        let mut coeff_dist : Vec<Z_q> = vec![];
        for i in 0..23 {
            coeff_dist.push(Z_q::from(0));
        }
        for i in 0..31 {
            coeff_dist.push(Z_q::from(1));
        }
        for i in 0..10 {
            coeff_dist.push(Z_q::from(2));
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

    pub fn generate_psi(&self) -> Vec<Z_q> {
        let mut rng = rand::thread_rng();
        let mut psi_k: Vec<Z_q> = Vec::new();
        for i in 0..L {
            psi_k.push(Z_q::from(rng.gen_range(0..Q)));
        }
        psi_k
    }
            

    pub fn generate_omega(&self) -> Vec<Z_q> {
        let mut rng = rand::thread_rng();
        let mut omega_k: Vec<Z_q> = Vec::new();
        for i in 0..256 {
            omega_k.push(Z_q::from(rng.gen_range(0..Q)));
        }
        //self.omega_k = Some(&omega_k);
        omega_k 
    }


    pub fn verify_b_prime_prime(&self, b_prime_prime_k : &R_q, omega_k: &Vec<Z_q>, psi_k : &Vec<Z_q>, projection: &Vec<i128>) -> () {
        // TODO again column vs row not sure.
        // Also self vs no self keyword not sure.
        let prod = Z_q::from(vec_inner_product(&Z_q::lift_inv(omega_k), projection));
        let mut sum = Z_q::zero();
        for i in 0..L {
            sum += &psi_k[i] * &(self.b_prime.as_ref().unwrap()[i]);
        }
        let check_val = prod + sum;

        // check that the constant term is equal to the above stuff.
        assert!(b_prime_prime_k.eval(Z_q::zero()) == check_val, "verify_b_prime_prime check failed: b_prime_prime_k constant coefficient is: {}, and check_val is: {}", b_prime_prime_k.eval(Z_q::zero()), check_val);
    }


    pub fn sample_jl_projection(&mut self) -> &Array2<i128> {
        let mut rng = rand::thread_rng();
        let choices = [-1,0,1];
        let weights = [0.25, 0.5, 0.25];
        let dist = WeightedIndex::new(&weights).unwrap();


        let mut Pi_i : Array2<i128> = Array2::zeros((256, N*(D as usize)));

        for ((i, j), value) in Pi_i.indexed_iter_mut() {
            *value = choices[dist.sample(&mut rng)] as i128;
        }

        self.Pi.as_mut().unwrap().push(Pi_i);

        // TODO don't entirely understand the functionality of this line.. but seems to work.
        let Pi_ref = self.Pi.as_ref().and_then(|pi| pi.last()).unwrap();
        Pi_ref
    }


    pub fn valid_projection(&mut self, projection: &Vec<i128>) -> bool {
        let val : f64 = 128.;
        let norm = l2_norm(projection);
        println!("TOTAL NORM OF JL PROJECTION: {}, {}", norm, val.sqrt()*(BETA_BOUND as f64));
        return norm <= val.sqrt()*(BETA_BOUND as f64);
    }

}

