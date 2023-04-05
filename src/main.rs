use polynomial::Polynomial;
extern crate nalgebra as na;
use na::{Matrix, DMatrix};
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;
use na::base::DVector;
use ndarray::{Array2, concatenate};

pub fn generate_polynomial() -> Polynomial<i32> {
    let mut poly = Polynomial::new();
    let num_terms = rng.gen_range(1, max_degree + 2);
    for degree in 0..num_terms {
        let coeff = rng.gen_range(0, q);
        p.add_coeff(degree, coeff);
    }
    poly
}


// Computes the ring inverse of a polynomial element of R_q
pub fn conjugation_automorphism(a : Polynomial) -> Polynomial {



}





// compute the dot product between two lists of 
pub fn polynomial_vec_inner_product(v1, v2) -> Polynomial<i32> {
    let result = Polynomial::new();
    for i in v1.length() {
        let product = v1[i].mul(&v2[i]);
        result.add(product);
    }
    result
}



pub fn create_proof() -> Result<Proof> {

    // compute witness: r vectors s_1, ..., s_r in R_q^n such that
    // f(s_1, .., s_r) = 0


    // constraints of the form
    // f(s) = sum a_ij<s_i, s_j> + sum <phi_i, s_i> + b = 0

    // a_ij, b \in R_q
    // R_q = Z_q[x]/(x^d+1) "polynomials with prime order integer coefficients mod x^d+1"
    // e.g., p = 7, x^3+1


    // TODO: these values are totally random at the moment

    // q: modulus of the ring of integers
    let q : i32 = 128;
    
    // constants: totally random for now and should be changed later
    let n : i32 = 128;
    let r : i32 = 20;


    // setting bound for SIS 2-norm
    let beta_bound : i32 = 65536

    let S= Array2::from_elem((n,r), Polynomial::new()); 

    // Random Generation of Witness matrix S
    
    // total summed norm:
    let mut norm_sum = 0;


    for i in 1..r {
        for j in 1..n {
            let s_ij = generate_polynomial();
            S[[i,j]] = s_ij;
            // add norm to norm sum
            norm_sum += i32::pow(s_ij.norm(), 2);
        }
    }

    // if too big, scale down:
    if norm_sum > i32::pow(beta_bound,2) {
        let scale_factor = i32::pow(beta_bound,2) / norm_sum;
        // TODO for all polynomials: *p *= scale;
    }



    // random generation of the symmetric polynomial matrix A
    let A = Array2::from_elem((r,r), Polynomial::new()); 
    for i in 0..r {
        for j in 0..r {
            let a_ij = generate_polynomial();
            if A[[i,j]] == Polynomial::new() {
                A[[i,j]] = a_ij;
                A[[j,i]] = a_ij;
            }
        }
    }

    // random generation of random polynomial matrix Phi
    let Phi = Array2::from_elem((n,r), Polynomial::new()); 
    for i in 0..n {
        for j in 0..r {
            let phi_ij = generate_polynomial();
            Phi[[i,j]] = phi_ij;
        }
    }


    // compute a_product and phi_product:
    let a_product : Polynomial = Polynomial::new();
    for i in 1..r {
        for j in 1..r {
            let inner_prod = polynomial_vec_inner_product(S.column_mut(i), S.column_mut(j));
            let prod = A[[i,j]].mul(&inner_prod);
            a_product.add(&prod);
        }
    }

    let phi_product : Polynomial = Polynomial::new();
    for i in 1..r {
        let inner_prod = polynomial_vec_inner_product(Phi.column_mut(i), S.column_mut(i));
        phi_product.add(&inner_prod);
    }

    let b : Polynomial = a_product.add(&phi_product);

    let t = Array2::<i32>::zeros((m,r)); // Inner commitment vector 


    // Compute inner Ajtai commitments
    // t_i = As_i \in R_q^\kappa (should just be m-dimensional)
    for i in 0..r {
        // TODO this only works if we switch to using MatrixMN instead of Array2.
        let t_i = A * s.column(i);
        t.column_mut(i).assign(&t_i);
    }

    // Compute outer commitment u_1
    let u_1 = b * t;

    // JL projection


    let projected_vec = jl_projection(&s[..]);


    // TODO: Fiat-Shamir version, full interactive or both?


}

// q: the value by which all integer coeffients are modded by
let q = 


// reduce from an R1CS instance (relatively easy to construct)
//  to R (module-SIS target relation, very difficult to construct)
pub fn binary_r1cs_reduction(A, B, C, w) {
    // First, check that this is in fact a valid R1CS instance:
    let a = A.dot(&w);
    let b = B.dot(&w);
    let c = C.dot(&w);
    assert!(((a.dot(&b) - c) == 0), "Invalid R1CS instance");

    // TODO later also check that this is binary


    // Sanity check: the verifier selects a polynomial A, and you're evaluating on the concatenation below
    let t = polynomial_a.evaluate(Array::concatenate[Axis(1), &[&a,&b,&c,&w]].unwrap()) % q;

    // Send to verifier. Verifier sends back 

    let phi = 

    for i in 1..l {
        phi_i = 

    }


    let a_inv = conjugation_automorphism(a);
    let b_inv = conjugation_automorphism(b);
    let c_inv = conjugation_automorphism(c);
    let w_inv = conjugation_automorphism(w);
}


// these have slightly different algorithms
pub fn arithmetic_r1cs_reduction() {



}




// generates and returns the SIS vector s
pub fn generate_witness() {
    






}




pub fn jl_projection(w: &[i32]) -> Vec<i32> {

    // proving knowledge of a long vector w in Z^d without revealing it

    // Let verifier sample a random linear map Pi: Z^d \to Z^256
    // entries of Pi are independent and equal to -1,0,1 with probabilities 1/4, 1/2, 1/4

    // sample random map Pi

    let d = w.len();


    let mut projection_vec = vec![0; 256];

    for i in 0..255 {
        let mut random_buffer = vec![0;d];
        for j in 0..(d-1) {
            if rand::random() {
                if rand::random() {
                    random_buffer[i] = 1;
                }
                else {
                    random_buffer[i] = -1;
                }
            }
            else {
                random_buffer[i] = 0;
            }
        }
        let v1 = DVector::from_vec(random_buffer);
        let v2 = DVector::from_row_slice(w);
        projection_vec[i] = v1.dot(&v2);
    }
    return projection_vec;
}





/*
pub fn verify_proof() -> Result<(), VerificationError> {
    // final check
    if true {
        Ok(())
    }

    else {
        Err(VerificationError::InvalidProof)
    }
}
*/


fn main() {

    // Test JL implementation with a random 1024-dimensional vector
    let mut rng = rand::rngs::StdRng::seed_from_u64(10); // set a seed for reproducibility
    let dist = Uniform::new(-100, 100);
    let my_vec: Vec<i32> = (0..1024).map(|_| rng.sample(&dist)).collect();
    let projected_vec = jl_projection(&my_vec[..]);

    let og_norm = f64::from(my_vec.iter().map(|x| x * x).sum::<i32>()).sqrt();
    let projected_norm = f64::from(projected_vec.iter().map(|x| x * x).sum::<i32>()).sqrt();

    let x: f64 = 30.0;

    println!("Original vector norm: {:.4}", og_norm);
    println!("Projected vector norm: {:.4}", projected_norm);
    println!("Projected scaled by sqrt(128): {:.4}", projected_norm*x.sqrt());
}
