pub mod util;
pub mod proofgen;
pub mod constants;
pub mod verification;
pub mod structs;


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
use crate::constants::*;
use crate::structs::*;
use crate::verification::*;

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

    // test norm computation of Array2 of polynomials
    let witness_norm : f64 = compute_total_norm(S);
    println!("Computed norm of witness matrix: {:.64}", witness_norm);

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
    println!("Welcome to the LaBRADOR Proof System!");
    println!("Generating Witness Matrix S");
    let S = generate_witness();

    println!("Generating Common Reference String (CRS)");
    let crs = CRS::new();

    println!("Generating State");
    let st = State::new(&S);

    let mut verifier = Verifier::new();
    let mut prover = Prover::new(&S, &mut verifier);

    println!("Generating proof..");
    let proof_transcript : Transcript = prover.proof_gen(&st, &crs);
    println!("Generated proof!");

    let res : bool = verifier.verify(&st, proof_transcript, &crs);
}

