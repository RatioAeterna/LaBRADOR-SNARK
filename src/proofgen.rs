use ndarray::{Array2, Ix2, concatenate, Axis};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;

pub fn proof_gen(S: Array2<Polynomial<i64>>) {
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
            let inner_prod = polynomial_vec_inner_product(S.column(i).to_vec(), S.column(j).to_vec());
            let prod = A[[i,j]].clone() * inner_prod;
            a_product = a_product.clone() + prod;
        }
    }

    let mut phi_product : Polynomial<i64> = Polynomial::new(vec![]);
    for i in 1..R {
        let inner_prod = polynomial_vec_inner_product(Phi.column(i).to_vec(), S.column(i).to_vec());
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
        let t_i = polynomial_matrix_product(A.clone(), S.t().column(i).to_owned().insert_axis(Axis(1))); // we need to compute the transpose of S for this product to work (for some reason) // we need to compute the transpose of S for this product to work (for some reason)
        println!("A dim {:?}", A.dim());
        println!("S.t() dim {:?}", S.t().dim());
        println!("S col dim {:?}", S.t().column(i).dim());
        println!("t_i dim: {:?}", t_i.dim());
        t.column_mut(i).assign(&t_i.remove_axis(Axis(1)));
    }
    println!("Computed Inner Ajtai Commitments!");
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
            norm_sum += poly_norm(s_ij);
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