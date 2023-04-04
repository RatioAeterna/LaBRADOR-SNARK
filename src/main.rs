use polynomial::Polynomial;
extern crate nalgebra as na;
use na::{Matrix, DMatrix};
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;
use na::base::DVector;
use ndarray::{Array2, concatenate};

pub fn create_proof() -> Result<Proof> {

    // compute witness: r vectors s_1, ..., s_r in R_q^n such that
    // f(s_1, .., s_r) = 0


    // constraints of the form
    // f(s) = sum a_ij<s_i, s_j> + sum <phi_i, s_i> + b = 0

    // a_ij, b \in R_q
    // R_q = Z_q[x]/(x^d+1) "polynomials with prime order integer coefficients mod x^d+1"
    // e.g., p = 7, x^3+1

    let b = Polynomial::new(vec)

    let s = Array2::<i32>::zeros((n,r)); // don't know what the witness is yet

    let t = Array2::<i32>::zeros((m,r)); // Inner commitment vector 


    // Compute inner Ajtai commitments
    // t_i = As_i \in R_q^\kappa (should just be m-dimensional)
    for i in 0..r {
        let t_i = A.matmul(&s.column(i)); 
        t.column_mut(i).assign(&t_i);
    }

    // Compute outer commitment u_1
    let u_1 = b.dot(&t);


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






}


// these have slightly different algorithms
pub fn arithmetic_r1cs_reduction() {



}





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
