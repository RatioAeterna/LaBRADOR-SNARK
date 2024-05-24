use criterion::{black_box, criterion_group, criterion_main, Criterion};
use labrador_snark::algebraic::*;
use labrador_snark::constants::*;
use labrador_snark::util::*;
//use std::sync::atomic::{AtomicBool, Ordering, AtomicOrdering};
use std::sync::atomic::Ordering as AtomicOrdering;

fn polynomial_multiplication(lhs: &Rq, rhs: &Rq) -> Rq {
    return lhs * rhs;
}

// NOTE: this function doesn't really account for overflow, so it's not really meant to be
// mathematically precise -- just for algorithm runtime estimation purposes
fn polynomial_multiplication_raw(lhs: &[u64], rhs: &[u64]) -> Vec<u64> {
    if NTT_ENABLED.load(AtomicOrdering::SeqCst) {
        // Perform NTT-based polynomial multiplication
        let mut prod: Vec<u64> = vec![0; D as usize];
        (*PLAN).negacyclic_polymul(&mut prod, &lhs, &rhs);
        prod
    } 
    else {
        // Perform classic polynomial multiplication
        let mut result = vec![0; lhs.len() + rhs.len() - 1];
        for (i, &a) in lhs.iter().enumerate() {
            for (j, &b) in rhs.iter().enumerate() {
                result[i + j] += a * b;
            }
        }
        // Next, add in reduction

        let mut reduced_result : Vec<u64> = vec![0; D];

        let data_len = result.len();
        for deg in (D as usize)..data_len {
            let term: u64 = (&result)[deg];
            // reduce the degree (which is > D) by dividing it by D
            let factor = (-1i128).pow((deg / D as usize) as u32); // TODO too big ints?
            let new_deg = deg % (D as usize);

            reduced_result[new_deg] = u64::from(Zq::from(term) * factor);
        }
        for deg in 0..(D as usize) {
            reduced_result[deg] = reduced_result[deg] + result[deg]; 
        }

        reduced_result
    }
}


fn bench_ntt_speed(c: &mut Criterion) {
    c.bench_function("poly_mul_classic", |b| {
        NTT_ENABLED.store(false, AtomicOrdering::SeqCst);
        let lhs: Rq = generate_polynomial(*Q, D);
        let rhs: Rq = generate_polynomial(*Q, D);
        b.iter(|| polynomial_multiplication(black_box(&lhs), black_box(&rhs)))
    });

    c.bench_function("poly_mul_ntt", |b| {
        NTT_ENABLED.store(true, AtomicOrdering::SeqCst);
        let lhs: Rq = generate_polynomial(*Q, D);
        let rhs: Rq = generate_polynomial(*Q, D);
        b.iter(|| polynomial_multiplication(black_box(&lhs), black_box(&rhs)))
    });
}

fn bench_ntt_raw(c: &mut Criterion) {
    c.bench_function("poly_mul_classic_raw", |b| {
        NTT_ENABLED.store(false, AtomicOrdering::SeqCst);
        let lhs: Vec<u64> = generate_polynomial(*Q, D).raw_coeffs();
        let rhs: Vec<u64> = generate_polynomial(*Q, D).raw_coeffs();
        b.iter(|| polynomial_multiplication_raw(black_box(&lhs), black_box(&rhs)))
    });

    c.bench_function("poly_mul_ntt_raw", |b| {
        NTT_ENABLED.store(true, AtomicOrdering::SeqCst);
        let lhs: Vec<u64> = generate_polynomial(*Q, D).raw_coeffs();
        let rhs: Vec<u64> = generate_polynomial(*Q, D).raw_coeffs();
        b.iter(|| polynomial_multiplication_raw(black_box(&lhs), black_box(&rhs)))
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(1000);
    targets = bench_ntt_speed, bench_ntt_raw
}
criterion_main!(benches);
