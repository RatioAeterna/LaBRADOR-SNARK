use ndarray::Array2;
use rand::{Rng};
use rand::distributions::{Distribution, WeightedIndex};

use num_traits::Zero;
use crate::util::*;
use crate::algebraic::*;
use crate::constants::*;
use crate::structs::*;

pub struct Verifier {
    b_prime : Option<Vec<Zq>>,
}

impl Verifier {

    pub fn new(b_prime_k: Vec<Zq>) -> Self {
        Verifier {
            b_prime: Some(b_prime_k),
        }
    }


    pub fn verify(&self, st: &State, proof : &Transcript, crs : &CRS) -> bool {

        // LINE 1, LINE 2
        // Enumeration of the state and transcript, respectively, so we don't include here.
        let upper_bound : usize = std::cmp::min(K, (128.0f64 / (*Q as f64).log2()).ceil() as usize);

        println!("starting line 3");
        // LINE 3
        // Computing "a_prime_prime" matrices for all k up to upper_bound, storing those in a vector

        let mut a_prime_prime : Vec<Array2<Rq>> = vec![];
        for k in 0..upper_bound {
            let mut a_prime_prime_mat = Array2::from_elem((R,R), Rq::new(vec![])); 
            for i in 0..R {
                for j in 0..R {
                    let mut sum : Rq = Rq::new(vec![]);
                    for l in 0..L {
                        let scaled : Rq = scale_polynomial(&st.a_prime_k[l][[i,j]], &proof.psi[k][l]);
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
        let mut phi_prime_prime_k : Vec<Vec<Vec<Rq>>> = vec![];
        for k in 0..upper_bound {
            //println!("k={}", k);
            // contains all phi_prime_prime_i for this particular k
            let mut phi_prime_prime : Vec<Vec<Rq>> = vec![];
            for i in 0..R {
                // ith polynomial vec for this particular k
                let mut sum : Vec<Rq> = vec![];
                for l in 0..L {
                    // TODO yes, we'll make this faster eventually
                    //println!("k: {}, l: {}, phi_prime_k len: {}, psi len: {}, psi_0 len: {}", k, l, &st.phi_prime_k.len(), proof.psi.len(), proof.psi[0].len());
                    let prod = scale_poly_vec(&st.phi_prime_k[l].column(i).to_vec(), &proof.psi[k][l]);
                    if sum.len() == 0 {
                        sum = vec![Rq::zero(); prod.len()];
                    }

                    sum = add_poly_vec(&sum, &prod);
                }
                for j in 0..256 {
                    let bolded_pi_poly_vec : Vec<Rq> = concat_coeff_reduction(&proof.Pi_i_all[i].row(j).to_vec());
                    let conj = sigma_inv_vec(&bolded_pi_poly_vec);
                    let prod = scale_poly_vec(&conj, &proof.omega[k][j]);
                    sum = add_poly_vec(&sum, &prod);
                }
                phi_prime_prime.push(sum);
            }
            phi_prime_prime_k.push(phi_prime_prime);
        }

        println!("starting line 5");
        // LINE 5
        // Forming a single canonical matrix of "a_ij" polynomials
        let mut Aij = Array2::from_elem((R,R), Rq::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                // generate single a_ij
                let mut a_ij : Rq = Rq::new(vec![]);
                for k in 0..K {
                    let prod : Rq = &proof.alpha[k] * &st.a_k[k][[i,j]];
                    a_ij = a_ij + prod;
                }
                for k in 0..upper_bound {
                    let prod : Rq = &proof.beta[k] * &a_prime_prime[k][[i,j]];
                    a_ij = a_ij + prod;
                }
                Aij[[i,j]] = a_ij;
            }
        }

        println!("starting line 6");
        // LINE 6 
        // Forming a single canonical set of vectors "phi_i" for i in {0, ..., R-1}
        let mut phi : Vec<Vec<Rq>> = vec![];
        for i in 0..R {
            let mut phi_i : Vec<Rq> = vec![]; 
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

        println!("starting line 7");
        // LINE 7
        // Forming a single canonical "b" polynomial
        let mut b : Rq = Rq::new(vec![]);
        for k in 0..K {
            let prod : Rq = &proof.alpha[k] * &st.b_k[k];
            b = b + prod;
        }
        for k in 0..upper_bound {
            let prod : Rq = &proof.beta[k] * &proof.b_prime_prime[k];
            b = b + prod;
        }
        
        println!("starting line 8");
        // CHECK 8
        // check that g_{ij} == g_{ji} i.e., matrix Gij is symmetric
        // TODO is it faster to only check some values? I really don't think this makes a
        // difference.
        for i in 0..R {
            for j in 0..R {
                if proof.Gij[[i,j]] != proof.Gij[[j,i]] {
                    return false;
                }
            }
        }

        println!("starting line 9");
        // CHECK 9
        // check that h_{ij} == h{ji} i.e., matrix Hij is symmetric
        for i in 0..R {
            for j in 0..R {
                if proof.Hij[[i,j]] != proof.Hij[[j,i]] {
                    //println!("At i,j: {}", proof.Hij[[i,j]]);
                    //println!("At j,i: {}", proof.Hij[[j,i]]);
                    return false;
                }
            }
        }

        println!("starting line 10");

        // LINE 10
        // Decompose vec z into z = z^(0) + z^(1)b
        let z_decompositions : Vec<Vec<Rq>> = decompose_polynomial_vec(&proof.z, *B, 2);

        println!("starting line 11");

        // LINE 11
        // Decompose vec t_i the same way
        let mut t_decompositions: Vec<Vec<Vec<Rq>>> = vec![];
        for i in 0..R {
            let t_i_decomposed : Vec<Vec<Rq>> = decompose_polynomial_vec(&proof.t_i_all[i], *B_1, *T_1);
            t_decompositions.push(t_i_decomposed);
        } 
        
        println!("starting line 12");

        // LINE 12
        // Decompose matrix elements g_ij
        let mut Gij_decompositions = Array2::from_elem((R,R), vec![]); 
        for i in 0..R {
            for j in 0..R {
                //println!("Gij term at [i={},j={}]... {}", i,j,&proof.Gij[[i,j]]);
                let dec_gij: Vec<Rq> = decompose_polynomial(&proof.Gij[[i,j]], *B_2, *T_2);
                Gij_decompositions[[i,j]] = dec_gij;
                //println!("decomposition.... {:?}", Gij_decompositions[[i,j]]);
            }
        }


        println!("starting line 13");

        // LINE 13
        // Decompose matrix elements h_ij
        let mut Hij_decompositions = Array2::from_elem((R,R), vec![]); 
        for i in 0..R {
            for j in 0..R {
                let dec_hij: Vec<Rq> = decompose_polynomial(&proof.Hij[[i,j]], *B_1, *T_1);
                Hij_decompositions[[i,j]] = dec_hij;
            }
        }


        println!("starting line 14");

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
                //println!("Gij term at [i={},j={}]... {}", i,j,&proof.Gij[[i,j]]);
                for k in 0..(*T_2 as usize) {
                    sum += poly_norm(&Gij_decompositions[[i,j]][k]); 
                    println!("Gij decomposition at k={}, is: {}", k,&Gij_decompositions[[i,j]][k]);
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
        println!("sum: {}, should be less than... {}", sum, (*BETA_PRIME).powi(2));

        if sum > (*BETA_PRIME).powi(2) { return false; }


        println!("starting line 15");

        // CHECK 15
        let mut lhs = polynomial_matrix_product(&crs.A_mat, &vec_to_column_array(&proof.z)).column(0).to_vec();
        let mut rhs : Vec<Rq> = vec![];
        for i in 0..R {
            let prod = poly_by_poly_vec(&proof.c[i], &proof.t_i_all[i]);
            if rhs.len() == 0 {
                rhs = gen_empty_poly_vec(prod.len());
            }
            rhs = add_poly_vec(&rhs, &prod);
        }

        if lhs != rhs { return false;}


        println!("starting line 16");

        // CHECK 16
        let lhs : Rq = polynomial_vec_inner_product(&proof.z, &proof.z);
        let mut rhs : Rq = Rq::new(vec![]);

        for i in 0..R {
            for j in 0..R {
                rhs = rhs + (&proof.Gij[[i,j]] * &proof.c[i] * &proof.c[j]);
                //println!("rhs in comp: {}", rhs);
            }
        }

        if lhs != rhs {
            return false;
        }

        println!("starting line 17");


        // CHECK 17
        let mut lhs : Rq = Rq::new(vec![]);
        //println!("{} {}", &proof.c.len(), &phi.len());
        for i in 0..R {
            lhs = lhs + polynomial_vec_inner_product(&phi[i], &proof.z) * &proof.c[i];
        }
        let mut rhs : Rq = Rq::new(vec![]);
        for i in 0..R {
            for j in 0..R {
                rhs = rhs + (&proof.Hij[[i,j]] * &proof.c[i] * &proof.c[j]);
            }
        }

        if lhs != rhs { return false;}


        println!("starting line 18");
    
        // CHECK 18
        let mut s1 : Rq = Rq::new(vec![]);
        let mut s2 : Rq = Rq::new(vec![]);
        for i in 0..R {
            for j in 0..R {
                let prod : Rq = &Aij[[i,j]] * &proof.Gij[[i,j]];
                s1 = s1 + prod;
            }
            s2 = s2 + &proof.Hij[[i,i]];
        }
        // check to make sure this is the zero polynomial
        if (s1 + s2 - b) != Rq::new(vec![]) { return false; }


        println!("starting line 19");

        // CHECK 19
        //let mut u_1_candidate : Vec<Rq> = vec![Rq::new(vec![]); KAPPA_1 as usize];
        let mut lhs = gen_empty_poly_vec(KAPPA_1 as usize);
        for i in 0..R {
            for k in 0..(*T_1 as usize) {
                let B_ik = crs.B_mat.get(&(i,k)).unwrap();
                let col_mat = vec_to_column_array(&t_decompositions[i][k]);
                let prod = polynomial_matrix_product(&B_ik, &col_mat).column(0).to_vec();
                lhs = add_poly_vec(&prod, &lhs);
            }
        }

        let mut rhs = gen_empty_poly_vec(KAPPA_2 as usize);
        for i in 0..R {
            for j in i..R {
                for k in 0..(*T_2 as usize) {
                    let C_ijk = crs.C_mat.get(&(i,j,k)).unwrap().column(0).to_vec();
                    let poly = &Gij_decompositions[[i,j]][k];
                    let prod = poly_by_poly_vec(poly, &C_ijk);
                    rhs = add_poly_vec(&prod, &rhs);
                }
            }
        }
        let u_1_candidate : Vec<Rq> = add_poly_vec(&lhs, &rhs);

        //println!("PROOF u1: {:?}, COMPUTED u1... {:?}", &proof.u_1, &u_1_candidate);
        //println!("lhs: {:?}, rhs: {:?}", &lhs, &rhs);
        if &proof.u_1 != &u_1_candidate { return false; }

        println!("starting line 20");

        // CHECK 20
        let mut u_2_candidate : Vec<Rq> = vec![Rq::new(vec![]); KAPPA_2 as usize];
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
        if &proof.u_2 != &u_2_candidate { return false; }
        true
    }

    // TODO do we want to STORE alpha, beta in the Verifier struct?
    pub fn fetch_alpha(&self) -> Vec<Rq> {
        let mut alpha = vec![]; 
        for i in 0..K {
            alpha.push(generate_polynomial(*Q,D));
        }
        alpha
    }

    pub fn fetch_beta(&self) -> Vec<Rq> {
        let mut beta = vec![]; 
        let upper_bound : usize = std::cmp::min(K, (128.0f64 / (*Q as f64).log2()).ceil() as usize);
        for i in 0..upper_bound {
            beta.push(generate_polynomial(*Q,D));
        }
        beta
    }


    // fetch a challenge polynomial from the challenge space \mathcal{C} satisfying a number of
    // criteria
    pub fn fetch_challenge(&self) -> Rq {
        // particular challenge coefficient distribution described on page 6.
        let mut coeff_dist : Vec<Zq> = vec![];
        if D == 64 {
            for i in 0..23 {
                coeff_dist.push(Zq::from(0));
            }
            for i in 0..31 {
                coeff_dist.push(Zq::from(1));
            }
            for i in 0..10 {
                coeff_dist.push(Zq::from(2));
            }
        }
        else {
            //coeff_dist.push(Zq::from(1));
            coeff_dist.push(Zq::from(0));
            //coeff_dist.push(Zq::from(0));
        }

        let candidate = generate_polynomial_picky(*Q,D as usize, coeff_dist.clone());
        // TODO... I don't think this norm should be squared, as it is.. which would perhaps give you 71^2.. definitely fix this if
        // needed.
        //assert!(poly_norm(&candidate) == TAU, "Incorrect l2 norm of challenge polynomial");
    
        while operator_norm(&candidate) > T {
            //assert!(poly_norm(&candidate) == TAU, "Incorrect l2 norm of challenge polynomial");
            let candidate = generate_polynomial_picky(*Q,D as usize, coeff_dist.clone());
        }
        candidate
    }

    pub fn generate_psi(&self) -> Vec<Zq> {
        let mut rng = rand::thread_rng();
        let mut psi_k: Vec<Zq> = Vec::new();
        for i in 0..L {
            //psi_k.push(Zq::from(rng.gen_range(0..*Q)));
            //println!("PSI K: {:?}", &psi_k);
            // TODO for debug purposes
            psi_k.push(Zq::from(0));
        }
        psi_k
    }
            

    pub fn generate_omega(&self) -> Vec<Zq> {
        let mut rng = rand::thread_rng();
        let mut omega_k: Vec<Zq> = Vec::new();
        for i in 0..256 {
            omega_k.push(Zq::from(rng.gen_range(0..*Q)));
            // TODO for debug purposes
            //omega_k.push(Zq::from(0));
        }
        //self.omega_k = Some(&omega_k);
        omega_k 
    }


    pub fn fetch_alleged_b_prime_prime_cc(&self, omega_k: &Vec<Zq>, psi_k: &Vec<Zq>, projection: &Vec<Zq>) -> Zq {

        //let prod = Zq::from(vec_inner_product(&Zq::lift_inv(omega_k), projection));
        let prod = Zq::from(vec_inner_product_Zq(omega_k, projection));
        let mut sum = Zq::zero();
        for i in 0..L {
            sum += &psi_k[i] * &(self.b_prime.as_ref().unwrap()[i]);
        }
        let alleged_b_prime_prime = prod + sum;


        alleged_b_prime_prime
    }




    pub fn verify_b_prime_prime(&self, b_prime_prime_k : &Rq, omega_k: &Vec<Zq>, psi_k : &Vec<Zq>, projection: &Vec<Zq>) -> () {
        // TODO again column vs row not sure.
        // Also self vs no self keyword not sure.
        //let prod = Zq::from(vec_inner_product(&Zq::lift_inv(omega_k), projection));
        let prod = Zq::from(vec_inner_product_Zq(omega_k, projection));
        let mut sum = Zq::zero();
        for i in 0..L {
            sum += &psi_k[i] * &(self.b_prime.as_ref().unwrap()[i]);
        }
        let check_val = prod + sum;

        // check that the constant term is equal to the above stuff.
        assert!(b_prime_prime_k.eval(Zq::zero()) == check_val, "verify_b_prime_prime check failed: b_prime_prime_k constant coefficient is: {}, and check_val is: {}", b_prime_prime_k.eval(Zq::zero()), check_val);
    }


    pub fn sample_jl_projection(&mut self) -> Array2<i128> {
        let mut rng = rand::thread_rng();
        let choices = [-1,0,1];
        let weights = [0.25, 0.5, 0.25];
        let dist = WeightedIndex::new(&weights).unwrap();


        let mut Pi_i : Array2<i128> = Array2::zeros((256, N*(D as usize)));

        for ((i, j), value) in Pi_i.indexed_iter_mut() {
            *value = choices[dist.sample(&mut rng)] as i128;
        }

        Pi_i
    }


    pub fn valid_projection(&mut self, projection: &Vec<i128>) -> bool {
        let val : f64 = 128.;
        let norm = l2_norm(projection);
        println!("TOTAL NORM OF JL PROJECTION: {}, {}", norm, val.sqrt()*(*BETA_BOUND as f64));
        return norm <= val.sqrt()*(*BETA_BOUND as f64);
    }

}

