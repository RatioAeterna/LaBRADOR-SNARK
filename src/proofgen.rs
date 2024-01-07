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

    pub fn proof_gen(&self, verifier : Verifier, crs : CRS) {
        //let t = Array2::<i64>::zeros((M,R)); // Inner commitment vector 
        //let t = Array2<Polynomial<i64>>::zeros((R,N)); // Inner commitment vector ( After looking at this some more, I think this is the actual dimension. Unclear, however. )
        let zero_poly = Polynomial::new(vec![0i64]);
        let mut t = Array2::from_elem((R,N), zero_poly);

        // Compute inner Ajtai commitments
        // t_i = As_i \in R_q^\kappa (should just be m-dimensional)
        for i in 0..R {
            //let t_i = A.clone().dot(S.t().column(i)); // we need to compute the transpose of S for this product to work (for some reason) // we need to compute the transpose of S for this product to work (for some reason)
            let t_i = polynomial_matrix_product(crs.A.clone(), self.S.t().column(i).to_owned().insert_axis(Axis(1))); // we need to compute the transpose of S for this product to work (for some reason) // we need to compute the transpose of S for this product to work (for some reason)
            println!("A dim {:?}", crs.A.dim());
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

                let res : Polynomial<i64> = polynomial_vec_inner_product(col_i, col_j);
                if Gij[[i,j]] == Polynomial::new(vec![]) {
                    Gij[[i,j]] = res;
                }
            }
        }

        let lhs = gen_empty_poly_vec(KAPPA_1);
        for i in 1..(R+1) {

            let t_i : Vec<Vec<Polynomial<i64>>> = decompose_polynomial_vec(t.column(i).to_vec(), B_1, T_1);
            for k in 0..(T_1) {
                let B_ik = crs.B.get((i,k)).unwrap();
                let t_mat = vec_to_column_array(t_i[k]);
                let prod = polynomial_matrix_product(B_ik, t_mat).to_vec();
                lhs = add_poly_vec(prod, lhs);
            }
        }

        let rhs = gen_empty_poly_vec(KAPPA_2);
        for i in 1..(R+1) {
            for j in i..(R+1) {
                let g_ij : Vec<Polynomial<i64>> = decompose_polynomial(Gij[[i,j]], B_2, T_2);
                for k in 0..(T_2) {
                    let C_ijk = crs.C.get((i,j,k)).unwrap().to_vec();
                    let g_ij_k : Polynomial<i64> = g_ij[k];
                    let prod = poly_by_poly_vec(g_ij_k, C_ijk);
                    rhs = add_poly_vec(prod, rhs);
                }
            }
        }

        // TODO investigate weird discrepancy between KAPPA_1 and KAPPA_2 dimension... 
        // FOR NOW, we assume they are equal rank. But I *think* you zero pad rhs if this isn't the
        // case.
        let u_1 : Vec<Polynomial<i64>> = add_poly_vec(lhs, rhs);











        // Next, compute JL projection
        let mut projection = jl_project(self.S.clone());
        while !verifier.valid_projection(projection) {
            println!("Verifier rejected projection!");
            projection = jl_project(self.S.clone());
        }
        println!("JL Projection complete and accepted");

        // First Aggregation step

        // TODO the paper mentions 'log' but it's not AT ALL clear whether this is base 2, 10, e, etc.
        // We will go with 10 for now.
        let upper_bound : usize = (128.0f64 / (q as f64).log10()).ceil();


        let psi : Vec<Vec<i64>> = vec![];
        let omega : Vec<Vec<i64>> = vec![];


        // TODO these comments aren't precise enough in dimensions
        // vector containing all phi'' for k in 1..(upper_bound), i in 0..(R-1)
        // TODO also perhaps switch this to an actual Vec<Array2<Polynomial<i64>>> later.
        let phi_prime_prime_vec : Vec<Vec<Vec<Polynomial<i64>>>> = vec![];

        // vector containing all b'' for k in 1..(upper_bound)
        let b_prime_prime_vec : Vec<Polynomial<i64>> = vec![];


        for k in 1..(upper_bound+1) {
            let psi_k = verifier.generate_psi_k();
            let omega_k = verifier.generate_omega_k();

            let phi_prime_prime_k : Vec<Vec<Polynomial<i64>>> = vec![];

            psi.push(psi_k);
            omega.push(omega_k);

            // NOTE: for our first basic rudimentary implementation, consider that
            // a_prime_ij is just going to be a_ij... given that our F' is just going to be F since
            // they both have zero constant coefficient. Simple F' for now.

            for i in 0..r {
                let phi_i_prime_prime : Vec<Polynoimal<i64>> = vec![];
                for j in 0..r {
                    let a_prime_ij : Polynomial<i64> = A[[i,j]];
                    let a_prime_prime_ij = multiply_poly_ints(a_prime_ij, psi_k);

                    // TODO is the Phi[i] indexing the right way in terms of column vs row? Just
                    // .col()
                    let phi_prime_i : Vec<Polynomial<i64>> = Phi[i].to_vec();


                    let mut lhs : Vec<Polynomial<i64>> = gen_empty_poly_vec(N);
                    for l in 1..(L+1) {
                        lhs = add_poly_vecs(multiply_poly_vec_ints(phi_prime_i , psi_k), lhs);
                    }

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

                    let phi_i_prime_prime = add_poly_vec_by_poly(lhs, rhs); 
                }
                phi_prime_prime_k.push(phi_i_prime_prime);
            }

            let mut lhs = Polynomial::new(vec![0i64]);
            for i in 0..R {
                for j in 0..R {
                    // TODO I think this is S column, but might be row. Double check later.
                    let prod = polynomial_vec_inner_product(S.column(i).to_vec(), S.column(j).to_vec());
                    let res = a_prime_prime[[i,j]] * prod;
                    lhs += res;
                }
            }

            let mut rhs = Polynomial::new(vec![0i64]);
            for i in 0..R {
                let res = polynomial_vec_inner_product(phi_prime_prime_k[i], S.column(i).to_vec());
                rhs += res;
            }
            let b_prime_prime_k = lhs + rhs;
            verifier.verify_b_prime_prime(b_prime_prime_k);

            phi_prime_prime_vec.push(phi_prime_prime_k);
            b_prime_prime_vec.push(b_prime_prime_k);
        }

        let alpha : Vec<Polynomial<i64>> = verifier.fetch_alpha();
        let beta : Vec<Polynomial<i64>> = verifier.fetch_beta();

        let phi_final : Vec<Vec<Polynomial<i64>> = vec![];
        for i in 0..R {
            let mut lhs : Vec<Polynomial<i64> = gen_empty_poly_vec(N);
            for k in 0..K {
                let res = poly_by_poly_vec(alpha[k], phi_k[k].column(i).to_vec());
                lhs = add_poly_vec(lhs, res);
            }
            let mut rhs : Vec<Polynomial<i64> = gen_empty_poly_vec(N);
            for k in 0..upper_bound {
                let res = poly_by_poly_vec(beta[k], phi_prime_prime_vec[k][i]);
                rhs = add_poly_vec(rhs , res);
            }
            phi_final.push(add_poly_vec(lhs, rhs));
        }


        // NOTE: This is *also* a symmetric matrix
        let mut Hij = Array2::from_elem((R,R), Polynomial::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let s_i = self.S.t().column(i).to_vec();
                let s_j = self.S.t().column(j).to_vec();

                let phi_i = phi_final[i];
                let phi_j = phi_final[j];

                let sum = polynomial_vec_inner_product(phi_i, s_j) + polynomial_vec_inner_product(phi_j, s_i);
                let res = scale_polynomial(sum, 0.5);

                if Hij[[i,j]] == Polynomial::new(vec![]) {
                    Hij[[i,j]] = res;
                }
            }
        }

        // Compute u_2
        let mut u_2 : Vec<Polynomial<i64>> = vec![];
        for i in 1..(R+1) {
            for j in i..(R+1) {
                for k in 0..T_1 {
                    let D_ijk_vec = crs.D.get((i,j,k)).unwrap().col(0).to_vec();
                    // TODO I don't think we can just add as such.
                    let dec_hij = decompose_polynomial(Hij[[i,j]]);
                    u_2 += poly_by_poly_vec(dec_hij[k], D_ijk_vec); 
                }
            }
        }

        let mut z : Vec<Polynomial<i64>> = vec![];
        let mut c_vec : Vec<Polynomial<i64>> = vec![];

        for i in 0..R {
            let c_i = verifier.fetch_challenge();
            c_vec.push(c_i);
            z.push(poly_by_poly_vec(c_i, S.column(i).to_vec()));
        }

        // TODO fill the proof transcript with all the relevant data and return
        let proof_transcript : Transcript = Transcript { projection, psi, omega, alpha, beta, u_2, c_vec, z, Gij, Hij };
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
