pub mod util;
pub mod proofgen;


use polynomial::Polynomial;
extern crate nalgebra as na;
use na::{Matrix, DMatrix};
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;
use na::base::DVector;
use ndarray::{Array2, concatenate};
use ndarray_linalg::norm;
use crate::util::*;
use crate::proofgen::*;

/*

// Computes the ring inverse of a polynomial element of R_q
pub fn conjugation_automorphism(a : Polynomial) -> Polynomial {
    // first, check that X^d+1 splits into two irreducible factors mod q.

    // we can do this using the Extended Euclidean algorithm:

    // for f in R_q, we find "a" such that (f*a congruent to 1) % x^d+1

    // we need Bezout coefficients of the following form:
    // f*a + (x^d+1)*b = 1, with gcd(f, x^d+1) = 1

    // step 1: x^d+1 = f*x + r






    a.mul(b)



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

    let t = Array2::<i64>::zeros((m,r)); // Inner commitment vector 


    // Compute inner Ajtai commitments
    // t_i = As_i \in R_q^\kappa (should just be m-dimensional)
    for i in 0..r {
        // TODO this only works if we switch to using MatrixMN instead of Array2.
        let t_i = A * S.column(i);
        t.column_mut(i).assign(&t_i);
    }

    // Compute outer commitment u_1
    let u_1 = b * t;

    // JL projection


    let projected_vec = jl_projection(&S[..]);


    // verifier checks that the ||projected_vec|| <= sqrt(128)*Beta
    valid_projection: bool = compute_norm(&DVector::from_vec(projected_vec)) <= sqrt(128)*BETA_BOUND);
    while !valid_projection {
        println!("Invalid JL Projection: Prover requests a new projection matrix Pi");
        let projected_vec = jl_projection(&s[..]);
        valid_projection: bool = compute_norm(&DVector::from_vec(projected_vec)) <= sqrt(128)*BETA_BOUND);
    }

    println!("Valid JL Projection")



    // First aggregation step


    // verifier must compute psi and omega (k-length random generated for k in 1..ceil(128/log q))




    let L = 1; // should be |F'|, which is... 1, for now?

    let psi: Vec<Vec<i64>> = Vec::new();
    let omega: Vec<Vec<i64>> = Vec::new();

    for k in 0..ceil(128/ (q as f64).log2() as u64) {
        // random draw from (Z_q)^L
        let psi_k = random_sample_Z_q(L, q);
        let omega_k = random_sample_Z_q(256, q);

        psi.push(psi_k);
        omega.push(omega_k);
    }








    // Amortization step

    // verifier sends c_i in C, subset of R_q. Prover computes the opening
    // z = c_1t_1 + ... + c_rt_r

    let z = Array2::<i64>::zeros((m,r)); // Random linear combination

    let C: Vec<Polynomial<i64>> = Vec::new(); // empty vector of challenge polynomials

    for i in 0..(r-1) {
        // multiply the challenge polynomial c_i by the VECTOR of polynomials t_i, sum all of these together.
        let c_i : Polynomial = generate_polynomial(); 
        let z_i = c_i * t.column(i); 
        z = z + z_i;
        C.push(c_i);
    }


    // Verifier checks the following:
    /*

    for i in 0..(r-1) {
        C.column(i)*t.column

    }


    assert!(A*z == );
    assert!(ndarray_linalg::norm::normalize(z) <= gamma);


    */


}


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

    /*

    let phi = 

    for i in 1..l {
        phi_i = 

    }

    */

    let a_inv = conjugation_automorphism(a);
    let b_inv = conjugation_automorphism(b);
    let c_inv = conjugation_automorphism(c);
    let w_inv = conjugation_automorphism(w);
}

// these have slightly different algorithms
pub fn arithmetic_r1cs_reduction() {



}






pub fn jl_projection(w: &[i64]) -> Vec<i64> {

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


*/

/*
TODO: better organization for testing later.
Right now we just test a bunch of util / misc. functions here haphazardly.
*/

#[test]
fn util_test() {
    // test polynomial generation
    let p1 = generate_polynomial(Q, D);

    println!("Random polynomial:");
    println!("{}", p1.pretty("x"));


    // test random sampling Z_q
    let sample = random_sample_Z_q(Q, 256);
    println!("String of integers mod q:");
    println!("{:?}", sample);

    // test scaling polynomial by a float < 1
    let scale: f32 = 0.5;
    let p2 = scale_polynomial(p1, scale);
    println!("Polynomial scaled by {:.32}:", scale);
    println!("{}", p2.pretty("x"));

    let S = generate_witness();
    println!("Generated witness matrix S:");
    for row in S.rows() {
        for poly in row {
            println!("{}", poly.pretty("x"));
        }
    }

    println!("Testing polynomial vec inner product:");
    let mut vec_a : Vec<Polynomial<i64>> = vec![];
    let mut vec_b : Vec<Polynomial<i64>> = vec![];
    for i in 0..5 {
        let mut new_poly = Polynomial::new(vec![0i64,1,2,3]);
        vec_a.push(new_poly.clone());
        vec_b.push(new_poly);
    }

    let product : Polynomial<i64> = Polynomial::new(vec![0i64,1,2,3]) * Polynomial::new(vec![0i64,1,2,3]) * Polynomial::new(vec![5i64]);

    let prod_poly = polynomial_vec_inner_product(vec_a, vec_b);
    assert!(prod_poly == product, "Polynomial inner products are broken.");
    println!("Polynomial products are working!");
    println!("{}", prod_poly.pretty("x"));
}

fn main() {

    println!("Hello, world!");
    let S = generate_witness();
    println!("Generated witness matrix S:");
    for row in S.rows() {
        for poly in row {
            println!("{}", poly.pretty("x"));
        }
    }
    println!("Generating proof..");
    proof_gen(S);
    println!("Generated proof!");
    /*
    // test polynomial generation
    let p1 = util::generate_polynomial(Q, D);

    println!("Random polynomial:");
    println!("{}", p1.pretty("x"));


    // test random sampling Z_q
    let sample = util::random_sample_Z_q(Q, 256);
    println!("String of integers mod q:");
    println!("{:?}", sample);
    */


}

