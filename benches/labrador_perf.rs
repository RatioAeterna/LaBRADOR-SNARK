use criterion::{black_box, BenchmarkId, criterion_group, criterion_main, Criterion};
use labrador_snark::algebraic::*;
use labrador_snark::constants::*;
use labrador_snark::structs::*;
use labrador_snark::verification::*;
use labrador_snark::proofgen::*;
use labrador_snark::util::*;
use core::time::Duration;
use std::sync::atomic::{AtomicBool, Ordering};

fn bench_labrador_perf(c: &mut Criterion) {
    let mut group = c.benchmark_group(
        "Perf for lattice size 2^n:",
    );
    group.sample_size(10);  // Sets a smaller sample size
    group.warm_up_time(Duration::from_secs(1));  // Shorter warm-up time
    group.measurement_time(Duration::from_secs(5));  // Shorter measurement time
    
    let mut n : usize = 1;
    let mut r : usize = 2;

    // TODO ideally we should be going from 2^2 to 2^20+, but this is just to test functionality
    for size_pow in 2..20 {
        if (size_pow % 2) == 0 {
            n = n * 2;
        }
        else {
            r = r * 2;
        }

        let constants = RuntimeConstants::new(n, r);

        let witness = generate_witness(&constants);
        let crs = CRS::new(&constants);
        let st = State::new(&witness, &constants);

        let mut verifier = Verifier::new(st.b_prime_k.clone(), &constants);
        let mut prover = Prover::new(&witness, &verifier, &constants);

        group.bench_with_input(BenchmarkId::from_parameter(size_pow), &size_pow, |b, &size_pow| {
            b.iter(|| {
                let proof_transcript: Transcript = prover.proof_gen(&st, &crs);
                let res: bool = verifier.verify(&st, &proof_transcript, &crs);
                assert!(res, "Error: proof verification failed!");
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_labrador_perf);
criterion_main!(benches);
