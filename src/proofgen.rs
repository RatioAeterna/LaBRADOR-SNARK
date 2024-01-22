use ndarray::{Array2, Ix2, concatenate, Axis}; use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;
use crate::verification::*;
use crate::structs::*;


pub struct Prover<'a> {

    // prover fields
    S: &'a Array2<Polynomial<i64>>,

    // reference to verifier
    // TODO does this even need to be mutable????
    verifier: &'a mut Verifier,

}

impl<'a> Prover<'a> {


    pub fn new(S: &'a Array2<Polynomial<i64>>, verifier: &'a mut Verifier) -> Self {
        Prover {
            S,
            verifier,
        }
    }

    pub fn proof_gen(&mut self, st : &State, crs : &CRS) -> Transcript {
        //let t = Array2::<i64>::zeros((M,R)); // Inner commitment vector 
        //let t = Array2<Polynomial<i64>>::zeros((R,N)); // Inner commitment vector ( After looking at this some more, I think this is the actual dimension. Unclear, however. )
        let zero_poly = Polynomial::new(vec![0i64]);
        let mut t = Array2::from_elem((R,N), zero_poly);

        // Compute inner Ajtai commitments
        // t_i = As_i \in R_q^\kappa (should just be m-dimensional)
        for i in 0..R {
            //let t_i = A.clone().dot(S.t().column(i)); // we need to compute the transpose of S for this product to work (for some reason) // we need to compute the transpose of S for this product to work (for some reason)
            let t_i = polynomial_matrix_product(&crs.A_mat, &self.S.t().column(i).to_owned().insert_axis(Axis(1))); // we need to compute the transpose of S for this product to work (for some reason) // we need to compute the transpose of S for this product to work (for some reason)
            println!("A dim {:?}", crs.A_mat.dim());
            println!("S.t() dim {:?}", self.S.t().dim());
            println!("S col dim {:?}", self.S.t().column(i).dim());
            println!("t_i dim: {:?}", t_i.dim());
            t.column_mut(i).assign(&t_i.remove_axis(Axis(1)));
        }
        println!("Computed Inner Ajtai Commitments!");


        // NOTE: This is a symmetric matrix
        let mut Gij = Array2::from_elem((R,R), Polynomial::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let col_i = self.S.t().column(i).to_vec();
                let col_j = self.S.t().column(j).to_vec();

                let res : Polynomial<i64> = polynomial_vec_inner_product(&col_i, &col_j);
                if Gij[[i,j]] == Polynomial::new(vec![]) {
                    Gij[[i,j]] = res;
                }
            }
        }

        let mut lhs = gen_empty_poly_vec(KAPPA_1 as usize);
        for i in 1..(R+1) {

            let t_i : Vec<Vec<Polynomial<i64>>> = decompose_polynomial_vec(&t.column(i).to_vec(), *B_1, *T_1);
            for k in 0..(*T_1 as usize) {
                let B_ik = crs.B_mat.get(&(i,k)).unwrap();
                let t_mat = vec_to_column_array(&t_i[k]);
                // NOTE this line might look somewhat confusing, but consider that the matrix
                // product is just a KAPPA_1 x 1 Array2 (which we turn into a vec by taking that
                // first column
                let prod = polynomial_matrix_product(B_ik, &t_mat).column(0).to_vec();
                lhs = add_poly_vec(&prod, &lhs);
            }
        }

        let mut rhs = gen_empty_poly_vec(KAPPA_2 as usize);
        for i in 1..(R+1) {
            for j in i..(R+1) {
                let g_ij : Vec<Polynomial<i64>> = decompose_polynomial(&Gij[[i,j]], *B_2, *T_2);
                for k in 0..(*T_2 as usize) {
                    let C_ijk = crs.C_mat.get(&(i,j,k)).unwrap().column(0).to_vec();
                    let g_ij_k : &Polynomial<i64> = &g_ij[k];
                    let prod = poly_by_poly_vec(g_ij_k, &C_ijk);
                    rhs = add_poly_vec(&prod, &rhs);
                }
            }
        }

        // TODO investigate weird discrepancy between KAPPA_1 and KAPPA_2 dimension... 
        // FOR NOW, we assume they are equal rank. But I *think* you zero pad rhs if this isn't the
        // case.
        let u_1 : Vec<Polynomial<i64>> = add_poly_vec(&lhs, &rhs);





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


        let mut psi : Vec<Vec<i64>> = vec![];
        let mut omega : Vec<Vec<i64>> = vec![];


        // TODO these comments aren't precise enough in dimensions
        // vector containing all phi'' for k in 1..(upper_bound), i in 0..(R-1)
        // TODO also perhaps switch this to an actual Vec<Array2<Polynomial<i64>>> later.
        let mut phi_prime_prime_vec : Vec<Vec<Vec<Polynomial<i64>>>> = vec![];

        // vector containing all b'' for k in 1..(upper_bound)
        let mut b_prime_prime_vec : Vec<Polynomial<i64>> = vec![];


        for k in 1..(upper_bound+1) {
            let psi_k = self.verifier.generate_psi();
            let omega_k = self.verifier.generate_omega();

            let mut phi_prime_prime_k : Vec<Vec<Polynomial<i64>>> = vec![];

            psi.push(psi_k);
            omega.push(omega_k);

            // NOTE: for our first basic rudimentary implementation, consider that
            // a_prime_ij is just going to be a_ij... given that our F' is just going to be F since
            // they both have zero constant coefficient. Simple F' for now.

            // TODO obviously we'll iterate the below look across all L in the end,
            // AND a_prime_prime will have a (k) superscript... so this will be a full K length vec eventually
            let mut a_prime_prime : Array2<Polynomial<i64>> = Array2::zeros((R,R));

            for i in 0..R {
                let phi_i_prime_prime : Vec<Polynomial<i64>> = vec![];
                for j in 0..R {
                    let a_prime_ij : &Polynomial<i64> = &st.a_prime_k[0][[i,j]];
                    let a_prime_prime_ij : Polynomial<i64> = multiply_poly_ints(&a_prime_ij, &psi.last().unwrap());
                    a_prime_prime[[i,j]] = a_prime_prime_ij;

                    // TODO is the phi_k[0][i] indexing the right way in terms of column vs row? Just
                    // .col()
                    // TODO we're taking the k=0 index for now.. disregarding the rest
                    let phi_prime_i : Vec<Polynomial<i64>> = st.phi_k[0].column(i).to_vec();


                    let mut lhs : Vec<Polynomial<i64>> = gen_empty_poly_vec(N);
                    for l in 1..(L+1) {                                 //NOTE: This is "psi_k"
                        lhs = add_poly_vec(&multiply_poly_vec_ints(&phi_prime_i , &psi.last().unwrap()), &lhs);
                    }

                    // side note: consider that pi_i^(j) is the jth row of Pi_i for j = 1, ..., 256.
                    let Pi_i = self.verifier.get_Pi_i(i);

                    let mut rhs : Polynomial<i64> = Polynomial::new(vec![]);

                    // TODO I think this should be 0..256 or something depending on how for loops
                    // work
                    /*
                    for j in 1..256 {
                        let conj = sigma_inv(Pi_i.row(j));
                        let omega_k_j = omega_k[j];

                        let res = scale_polynomial(conj, omega_k_j as f32);
                        rhs = rhs + res;
                    }
                    */

                    let phi_i_prime_prime = add_poly_vec_by_poly(&lhs, &rhs); 
                }
                phi_prime_prime_k.push(phi_i_prime_prime);
            }

            let mut lhs = Polynomial::new(vec![0i64]);
            for i in 0..R {
                for j in 0..R {
                    // TODO I think this is S column, but might be row. Double check later.
                    let prod = polynomial_vec_inner_product(&self.S.column(i).to_vec(), &self.S.column(j).to_vec());
                    let res = &a_prime_prime[[i,j]] * prod;
                    lhs = lhs + res;
                }
            }

            let mut rhs = Polynomial::new(vec![0i64]);
            for i in 0..R {
                let res = polynomial_vec_inner_product(&phi_prime_prime_k[i], &self.S.column(i).to_vec());
                rhs = rhs + res;
            }
            let b_prime_prime_k = lhs + rhs;
            //self.verifier.verify_b_prime_prime(b_prime_prime_k);

            phi_prime_prime_vec.push(phi_prime_prime_k);
            b_prime_prime_vec.push(b_prime_prime_k);
        }

        let alpha : Vec<Polynomial<i64>> = self.verifier.fetch_alpha();
        let beta : Vec<Polynomial<i64>> = self.verifier.fetch_beta();

        let mut phi_final : Vec<Vec<Polynomial<i64>>> = vec![];
        for i in 0..R {
            let mut lhs : Vec<Polynomial<i64>> = gen_empty_poly_vec(N);
            for k in 0..K {
                let res = poly_by_poly_vec(&alpha[k], &st.phi_k[k].column(i).to_vec());
                lhs = add_poly_vec(&lhs, &res);
            }
            let mut rhs : Vec<Polynomial<i64>> = gen_empty_poly_vec(N);
            for k in 0..upper_bound {
                let res = poly_by_poly_vec(&beta[k], &phi_prime_prime_vec[k][i]);
                rhs = add_poly_vec(&rhs , &res);
            }
            phi_final.push(add_poly_vec(&lhs, &rhs));
        }


        // NOTE: This is *also* a symmetric matrix
        let mut Hij = Array2::from_elem((R,R), Polynomial::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let s_i = self.S.t().column(i).to_vec();
                let s_j = self.S.t().column(j).to_vec();

                let phi_i = &phi_final[i];
                let phi_j = &phi_final[j];

                let sum = polynomial_vec_inner_product(phi_i, &s_j) + polynomial_vec_inner_product(phi_j, &s_i);
                let res = scale_polynomial(&sum, 0.5);

                if Hij[[i,j]] == Polynomial::new(vec![]) {
                    Hij[[i,j]] = res;
                }
            }
        }

        // Compute u_2
        let mut u_2 : Vec<Polynomial<i64>> = vec![Polynomial::<i64>::new(vec![]); KAPPA_2 as usize];
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

        let mut z : Vec<Polynomial<i64>> = vec![Polynomial::new(vec![]); R as usize];
        let mut c_vec : Vec<Polynomial<i64>> = vec![];

        for i in 0..R {
            let c_i = self.verifier.fetch_challenge();
            c_vec.push(c_i);
            let prod : Vec<Polynomial<i64>> = poly_by_poly_vec(c_vec.last().as_ref().unwrap(), &self.S.column(i).to_vec());
            for j in 0..R {
                z[j] = &z[j] + &prod[j];
            }
        }

        // TODO fill the proof transcript with all the relevant data and return
        let proof_transcript : Transcript = Transcript { projection, psi, omega, alpha, beta, u_2, c : c_vec, z, Gij, Hij };
        proof_transcript
    }

    pub fn jl_project(&mut self) -> Vec<i64> {
        // verifier sends random matrices in {-1,0,1}
        let mut projection : Vec<i64> = vec![0 ; 256];
        for i in 0..R {
            let Pi_i = self.verifier.sample_jl_projection();
            let s_i_coeffs : Array2<i64> = vec_to_column_array(&witness_coeff_concat(&self.S.column(i).to_vec()));
            // NOTE: for reference, this is a 256x(ND) multiplied by an (ND)x1, giving a 256x1
            // which we turn into a vec
            let product = matmul(Pi_i, &s_i_coeffs).column(0).to_vec();
            projection = add_vecs(&projection, &product);
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
            norm_sum += poly_norm(&s_ij) as i64;
        }
    }

    // if too big, scale down:
    if norm_sum > (i64::pow(BETA_BOUND,2)) {
        let scale_factor : f32 = ((i64::pow(BETA_BOUND,2)) / norm_sum) as f32;
        // scale each polynomial in the matrix by scale factor
        for i in 0..N {
            for j in 0..R {
                S[[i,j]] = scale_polynomial(&S[[i,j]], scale_factor);
            }
        }
    }
    S
}
