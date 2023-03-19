use polynomial::Polynomial;
extern crate nalgebra as na;
use na::{Matrix, DMatrix};
use rand::prelude::*;

pub fn create_proof() -> Result<Proof> {

    // compute witness: r vectors s_1, ..., s_r in R_q^n such that
    // f(s_1, .., s_r) = 0


    // constraints of the form
    // f(s) = sum a_ij<s_i, s_j> + sum <phi_i, s_i> + b = 0

    // a_ij, b \in R_q
    // R_q = Z_q[x]/(x^d+1) "polynomials with prime order integer coefficients mod x^d+1"
    // e.g., p = 7, x^3+1



    let b = Polynomial::new(vec)
}


pub fn JL_projection() {

    // proving knowledge of a long vector w in Z^d without revealing it

    // Let verifier sample a random linear map Pi: Z^d \to Z^256
    // entries of Pi are independent and equal to -1,0,1 with probabilities 1/4, 1/2, 1/4

    // sample random map Pi

    w = vec![0;d];
    random_buffer = vec![0; d];
    for i in 0..d-1 {
        if rand::random() {
            if rand::random() {
                random_buffer[i] = 1;
            }
            else random_buffer[i] = -1;
        }
        else random_buffer[i] = 0;
    }
    let Pi = new Matrix<i64, 256, d, random_buffer>;
    return Pi * w;
}





pub fn verify_proof() -> Result<(), VerificationError> {
    // final check
    if true {
        Ok(())
    }

    else {
        Err(VerificationError::InvalidProof)
    }
}


fn main() {
    println!("Test");
}
