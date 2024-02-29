extern crate nalgebra as na;
use na::{Matrix, DMatrix};
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;
use na::base::DVector;
use ndarray::{Array2, concatenate};
use ndarray_linalg::norm;
use num_traits::One;
use num_traits::Zero;
use labrador_snark::algebraic::*;
use labrador_snark::util::*;
use labrador_snark::proofgen::*;
use labrador_snark::constants::*; 
use labrador_snark::structs::*;
use labrador_snark::verification::*;
use std::sync::atomic::{AtomicBool, Ordering};

/*
    TODO: better organization for testing later.
    Right now we just test a bunch of util / misc. functions here haphazardly.
*/
// TODO also look into using rustfmt with the LSP server, etc.

fn print_constants() {
    println!("Printing runtime-computed constants:\n");
    println!("STD: {}", *STD);
    println!("B: {}", *B);
    println!("B_1: {}", *B_1);
    println!("B_2: {}", *B_2);
    println!("T_1: {}", *T_1);
    println!("T_2: {}", *T_2);
    println!("GAMMA: {}", *GAMMA);
    println!("GAMMA_1: {}", *GAMMA_1);
    println!("GAMMA_2: {}", *GAMMA_2);
    println!("BETA_PRIME: {}", *BETA_PRIME);
}

fn main() {

    NTT_ENABLED.store(false, Ordering::SeqCst);

    println!("Welcome to the LaBRADOR Proof System!");

    print_constants();

    println!("Generating Witness Matrix S");
    let S = generate_witness();
    //println!("sanity check of witness vals: {}", S[[0,0]]);
    println!("sanity check of witness vals: {}", &S[[0,0]]);

    println!("Generating Common Reference String (CRS)");
    let crs = CRS::new();

    println!("Generating State");
    let st = State::new(&S);

    let mut verifier = Verifier::new(st.b_prime_k.clone());
    let mut prover = Prover::new(&S, &mut verifier);

    println!("Generating proof..");
    let proof_transcript : Transcript = prover.proof_gen(&st, &crs);
    println!("Generated proof!");

    println!("Verifying proof..");
    let res : bool = verifier.verify(&st, proof_transcript, &crs);
    assert!( res, "Error: Proof Verification Failed");
    println!("Success: Proof Verified!");
}

