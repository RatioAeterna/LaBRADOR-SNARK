use labrador_snark::proofgen::*;
use labrador_snark::constants::*; 
use labrador_snark::structs::*;
use labrador_snark::verification::*;
use std::sync::atomic::Ordering;
use clap::Parser;

// TODO also look into using rustfmt with the LSP server, etc.
fn print_constants() {
    println!("Printing runtime-computed constants:");
    println!("Q: {}", *Q);
    println!("BETA: {}", *BETA_BOUND);
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
    let verbose : bool = args.verbose;
    let ntt : bool = args.ntt;

    if verbose {
        let _ = VERBOSE.set(true);
    }
    else {
        let _ = VERBOSE.set(false);
    }

    if ntt {
        NTT_ENABLED.store(true, Ordering::SeqCst);
    }
    else {
        NTT_ENABLED.store(false, Ordering::SeqCst);
    }

    if is_verbose() {
        println!("Welcome to the LaBRADOR Proof System!");
        println!("=====================================\n");
        print_constants();
        println!("Generating Witness Matrix");
    }

    let witness = generate_witness();

    if is_verbose() {
        println!("sanity check of witness vals: {}", &witness[[0,0]]);
        println!("Generating Common Reference String (CRS)");
    }

    let crs = CRS::new();

    if is_verbose() { println!("Generating State"); }

    let st = State::new(&witness);

    let mut verifier = Verifier::new(st.b_prime_k.clone());
    let mut prover = Prover::new(&witness, &mut verifier);


    if is_verbose() { println!("Generating proof.."); }

    let proof_transcript : Transcript = prover.proof_gen(&st, &crs);

    if is_verbose() { println!("Generated proof!"); }

    if is_verbose() { println!("Verifying proof.."); }
    let res : bool = verifier.verify(&st, &proof_transcript, &crs);
    assert!( res, "Error: Proof Verification Failed");
    if is_verbose() { 
        println!("Success: Proof Verified!"); 
        println!("=========================");
        println!("Size of proof: {} KB", (proof_transcript.size_in_bytes() as f64)/ 1024.0);
    }
}

