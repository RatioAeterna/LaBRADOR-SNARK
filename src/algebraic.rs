use std::ops::Mul;
use ndarray::{Array2, Ix2, concatenate, Axis}; 
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use crate::constants::*;
use num_traits::Zero;

/*
 * Crate for various algebraic structure implementations
 *
 */

// The ring of integers modulo q (Q in constants.rs)
#[derive(Clone, Debug)]
struct Z_q {
    value: i128,
}

impl Z_q {
    fn new(value: i128) -> Self {
        Z_q { value : value % Q }
    }

}

impl Zero for Z_q {

    fn zero() -> Self {
        Z_q { value: 0 }
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }

}

impl std::ops::Add for Z_q {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::new(self.value + other.value)
    }
}

impl std::ops::Sub for Z_q {
    type Output = Self;

    // TODO is this correct? Double check
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.value - other.value)
    }
}

impl std::ops::Mul for Z_q {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        Self::new(self.value * other.value)
    }
}
// TODO do we need Div also? Probably not.

// Wrapper for Polynomial<Z_q> which specifically uses
// properties of R_q, aka Z[q]/(X^d+1)
struct R_q(Polynomial<Z_q>);

impl R_q {

    // TODO maybe overload so you can pass in refs instead..
    fn new(coefficients: Vec<Z_q>) -> Self {
        R_q(Polynomial::new(coefficients))
    }

}

impl std::ops::Add for &R_q {
    type Output = R_q;

    fn add(self, rhs: Self) -> Self::Output {
        let sum : Polynomial<Z_q> = &self.0 + &rhs.0;
        let sum_rq = R_q(sum);
        sum_rq
    }
}

impl std::ops::Sub for &R_q {
    type Output = R_q;

    fn sub(self, rhs: Self) -> Self::Output {
        let minus : Polynomial<Z_q> = &self.0 - &rhs.0;
        let minus_rq = R_q(minus);
        minus_rq
    }
}

impl std::ops::Mul for &R_q {

    type Output = R_q;
        
    fn mul(self, rhs: Self) -> Self::Output {
        // Custom multiplication logic.
        // We need to reduce by (X^d+1)

        // TODO add NTT logic here later.
        let prod : Polynomial<Z_q> = &self.0 * &rhs.0;
        let data_len = (&prod).data().len()+1;

        let mut prod_rq = R_q(prod);

        // REDUCE all terms of deg > 64
        for deg in (D as usize)+1..data_len {
            let term : Z_q = (&prod_rq).0.data()[deg].clone();
            let mut term_poly_data = vec![Z_q::zero(); deg+1];
            term_poly_data[deg] = term;
            let term_poly : R_q = R_q::new(term_poly_data);
            prod_rq = &prod_rq + &term_poly;
        }
        prod_rq
    }
}

