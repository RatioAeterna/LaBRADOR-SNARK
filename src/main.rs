use clap::Parser;
use labrador_snark::constants::*;
use labrador_snark::proofgen::*;
use labrador_snark::structs::*;
use labrador_snark::verification::*;
use labrador_snark::util::*;
use std::sync::atomic::Ordering;

// TODO also look into using rustfmt with the LSP server, etc.
fn print_constants(constants: &RuntimeConstants) {
    println!("Printing runtime-computed constants:");
    println!("Q: {}", *Q);
    println!("BETA: {}", constants.BETA_BOUND);
    println!("STD: {}", constants.STD);
    println!("B: {}", constants.B);
    println!("B_1: {}", constants.B_1);
    println!("B_2: {}", constants.B_2);
    println!("T_1: {}", constants.T_1);
    println!("T_2: {}", constants.T_2);
    println!("GAMMA: {}", constants.GAMMA);
    println!("GAMMA_1: {}", constants.GAMMA_1);
    println!("GAMMA_2: {}", constants.GAMMA_2);
    println!("BETA_PRIME: {}", constants.BETA_PRIME);
}

/// Robust implementation of the LaBRADOR Cryptographic Proof System
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Verbosity setting (default value of false / no printed output)
    #[arg(short, long)]
    verbose: bool,

    /// Enable Number Theoretic Transform for faster polynomial multiplication
    #[arg(short, long)]
    ntt: bool,
}

fn main() {
    let args = Args::parse();
    let verbose: bool = args.verbose;
    let ntt: bool = args.ntt;

    if verbose {
        let _ = VERBOSE.set(true);
    } else {
        let _ = VERBOSE.set(false);
    }

    if ntt {
        NTT_ENABLED.store(true, Ordering::SeqCst);
    } else {
        NTT_ENABLED.store(false, Ordering::SeqCst);
    }

    let n : usize = 2;
    let r : usize = 2;

    let constants = RuntimeConstants::new(n, r);

    if is_verbose() {
        println!("Welcome to the LaBRADOR Proof System!");
        println!("=====================================\n");
        print_constants(&constants);
        println!("Generating Witness Matrix");
    }

    let witness = generate_witness(&constants);


    if is_verbose() {
        println!("sanity check of witness vals: {}", &witness[[0, 0]]);
        println!("Generating Common Reference String (CRS)");
    }

    let crs = CRS::new(&constants);

    if is_verbose() {
        println!("Generating State");
    }

    let st = State::new(&witness, &constants);

    let mut verifier = Verifier::new(st.b_prime_k.clone(), &constants);
    let mut prover = Prover::new(&witness, &verifier, &constants);

    if is_verbose() {
        println!("Generating proof..");
    }

    let proof_transcript: Transcript = prover.proof_gen(&st, &crs);

    if is_verbose() {
        println!("Generated proof!");
    }

    if is_verbose() {
        println!("Verifying proof..");
    }
    let res: bool = verifier.verify(&st, &proof_transcript, &crs);
    assert!(res, "Error: Proof Verification Failed");
    if is_verbose() {
        println!("Success: Proof Verified!");
        println!("=========================");
        println!(
            "Size of proof: {} KB",
            (proof_transcript.size_in_bytes() as f64) / 1024.0
        );
    }
}
