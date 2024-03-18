use criterion::{black_box, criterion_group, criterion_main, Criterion};
use labrador_snark::constants::*;
use labrador_snark::algebraic::*;
use labrador_snark::util::*;
use std::sync::atomic::{AtomicBool, Ordering};

fn polynomial_multiplication(lhs : &Rq, rhs: &Rq) -> Rq {
    return lhs * rhs;
}

fn bench_ntt_speed(c: &mut Criterion) {
    c.bench_function("poly_mul_classic", |b| {
        NTT_ENABLED.store(false, Ordering::SeqCst);
        let lhs : Rq = generate_polynomial(*Q, D);
        let rhs : Rq = generate_polynomial(*Q, D);
        b.iter(|| polynomial_multiplication(black_box(&lhs), black_box(&rhs)))
    });

    c.bench_function("poly_mul_ntt", |b| {
        NTT_ENABLED.store(true, Ordering::SeqCst);
        let lhs : Rq = generate_polynomial(*Q, D);
        let rhs : Rq = generate_polynomial(*Q, D);
        b.iter(|| polynomial_multiplication(black_box(&lhs), black_box(&rhs)))
    });
}

criterion_group!(benches, bench_ntt_speed);
criterion_main!(benches);

