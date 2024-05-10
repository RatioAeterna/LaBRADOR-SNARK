use criterion::{black_box, criterion_group, criterion_main, Criterion};
use labrador_snark::algebraic::*;
use labrador_snark::constants::*;
use labrador_snark::util::*;
//use std::sync::atomic::{AtomicBool, Ordering, AtomicOrdering};
use std::sync::atomic::Ordering as AtomicOrdering;

fn polynomial_multiplication(lhs: &Rq, rhs: &Rq) -> Rq {
    return lhs * rhs;
}

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
        result
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

criterion_group!(benches, bench_ntt_speed, bench_ntt_raw);
criterion_main!(benches);
