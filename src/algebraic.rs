use std::ops::Mul;
use std::fmt;
use ndarray::{Array2, Ix2, concatenate, Axis}; 
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use crate::constants::*;
use num_traits::Zero;
use num_traits::One;
use std::cmp::Ordering;

/*
 * Crate for various algebraic structure implementations
 *
 */

// The ring of integers modulo q (Q in constants.rs)
#[derive(Clone, Copy, Debug, Eq, Ord, PartialOrd)]
pub struct Z_q {
    value: i128,
}

impl Z_q {
    pub fn new(value: i128) -> Self {
        Z_q { value : value % Q }
    }
    pub fn lift(vec : &Vec<i128>) -> Vec<Z_q> {
        let mut res_vec : Vec<Z_q> = vec![];

        for i in 0..vec.len() {
            res_vec.push(Z_q::from(vec[i]));
        }
        res_vec
    }
    pub fn lift_inv(vec : &Vec<Z_q>) -> Vec<i128> {
        let mut res_vec : Vec<i128> = vec![];
        for i in 0..vec.len() {
            res_vec.push(vec[i].value);
        }
        res_vec
    }

}

impl std::ops::Neg for Z_q {
    type Output = Z_q;
    fn neg(self) -> Z_q {
        Z_q {
            value: if self.value == 0 { 0 } else { Q - self.value },
        }
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

impl std::ops::Rem<i128> for Z_q {
    type Output = Z_q;

    fn rem(self, rhs: i128) -> Self::Output {
        Z_q::new(self.value % rhs)
    }
}






impl One for Z_q {
    fn one() -> Self {
        Z_q { value: 1 }
    }
}

impl std::iter::Sum for Z_q {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}






impl std::ops::Add for Z_q {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::new(self.value + other.value)
    }
}

impl std::ops::Add<i128> for Z_q {
    type Output = Self;

    fn add(self, other: i128) -> Self::Output {
        Self::new(self.value + other)
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

impl std::ops::Mul<i128> for Z_q {
    type Output = Self;

    fn mul(self, other: i128) -> Self::Output {
        Self::new(self.value * other)
    }
}

impl std::ops::Div<i128> for Z_q {
    type Output = Self;

    fn div(self, rhs: i128) -> Self::Output {
        Self::new(self.value / rhs)
    }
}


impl std::ops::Add for &Z_q {
    type Output = Z_q;

    fn add(self, other: Self) -> Self::Output {
        Z_q::new(self.value + other.value)
    }
}

impl std::ops::Sub for &Z_q {
    type Output = Z_q;

    // TODO is this correct? Double check
    fn sub(self, other: Self) -> Self::Output {
        Z_q::new(self.value - other.value)
    }
}

impl std::ops::Mul for &Z_q {
    type Output = Z_q;

    fn mul(self, other: Self) -> Self::Output {
        Z_q::new(self.value * other.value)
    }
}

impl std::ops::AddAssign<Z_q> for Z_q {
    fn add_assign(&mut self, other: Z_q) {
        self.value = (self.value + other.value) % Q;
    }
}

impl std::ops::SubAssign<Z_q> for Z_q {
    fn sub_assign(&mut self, other: Z_q) {
        self.value = (self.value - other.value) % Q;
    }
}



impl fmt::Display for Z_q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl From<&Z_q> for usize {
    fn from(item: &Z_q) -> Self {
        item.value as usize
    }
}


impl From<Z_q> for f32 {
    fn from(item: Z_q) -> Self {
        item.value as f32
    }
}

impl From<Z_q> for f64 {
    fn from(item: Z_q) -> Self {
        item.value as f64
    }
}

impl From<f32> for Z_q {
    fn from(item: f32) -> Self {
        Self::new(item as i128)
    }
}


impl From<f64> for Z_q {
    fn from(item: f64) -> Self {
        Self::new(item as i128)
    }
}

impl From<i32> for Z_q {
    fn from(item: i32) -> Self {
        Self::new(item as i128)
    }
}

impl From<i128> for Z_q {
    fn from(item: i128) -> Self {
        Self::new(item)
    }
}

impl PartialEq for Z_q {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl PartialEq<i128> for Z_q {
    fn eq(&self, other: &i128) -> bool {
        self.value == (other % Q)
    }
}

impl PartialOrd<i128> for Z_q {
    fn partial_cmp(&self, other: &i128) -> Option<Ordering> {
        self.value.partial_cmp(&(other % Q))
    }
}



// Wrapper for Polynomial<Z_q> which specifically uses
// properties of R_q, aka Z[q]/(X^d+1)
#[derive(Clone, Debug)]
pub struct R_q(Polynomial<Z_q>);

impl R_q {
    // TODO maybe overload so you can pass in refs instead..
    pub fn new(coefficients: Vec<Z_q>) -> Self {
        R_q(Polynomial::new(coefficients))
    }

    pub fn data_vec(&self) -> Vec<Z_q> {
        self.0.data().to_vec()
    }

    pub fn eval(&self,  x: Z_q) -> Z_q {
        self.0.eval(x)
    }


}



impl Zero for R_q {

    fn zero() -> Self {
        R_q::new(vec![])
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

}

impl std::ops::Add for R_q {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let sum : Polynomial<Z_q> = self.0 + rhs.0;
        let sum_rq = R_q(sum);
        sum_rq
    }
}

impl std::ops::Add<&R_q> for R_q {
    type Output = Self;

    fn add(self, rhs: &R_q) -> Self::Output {
        let sum : Polynomial<Z_q> = self.0 + &rhs.0;
        let sum_rq = R_q(sum);
        sum_rq
    }
}



impl std::ops::Sub for R_q {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let sum : Polynomial<Z_q> = self.0 - rhs.0;
        let sum_rq = R_q(sum);
        sum_rq
    }
}

impl std::ops::Mul for R_q {
    type Output = Self;

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
            // reduce the degree (which is > 64) by dividing it by 64
            let factor = (-1i128).pow((deg / 64) as u32); // TODO too big ints?
            let new_deg = deg % 64;

            let mut term_poly_data = vec![Z_q::zero(); deg+1];
            term_poly_data[new_deg] = term*factor;
            let term_poly : R_q = R_q::new(term_poly_data);
            prod_rq = &prod_rq + &term_poly;
        }
        prod_rq
    }
}

impl std::ops::Mul<&R_q> for R_q {
    type Output = Self;

    fn mul(self, rhs: &R_q) -> Self::Output {
        // Custom multiplication logic.
        // We need to reduce by (X^d+1)

        // TODO add NTT logic here later.
        let prod : Polynomial<Z_q> = &self.0 * &rhs.0;
        let data_len = (&prod).data().len()+1;

        let mut prod_rq = R_q(prod);

        // REDUCE all terms of deg > 64
        for deg in (D as usize)+1..data_len {
            let term : Z_q = (&prod_rq).0.data()[deg].clone();
            // reduce the degree (which is > 64) by dividing it by 64
            let factor = (-1i128).pow((deg / 64) as u32); // TODO too big ints?
            let new_deg = deg % 64;

            let mut term_poly_data = vec![Z_q::zero(); deg+1];
            term_poly_data[new_deg] = term*factor;
            let term_poly : R_q = R_q::new(term_poly_data);
            prod_rq = &prod_rq + &term_poly;
        }
        prod_rq
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

impl std::ops::Add<R_q> for &R_q {
    type Output = R_q;

    fn add(self, rhs: R_q) -> Self::Output {
        let sum : Polynomial<Z_q> = &self.0 + rhs.0;
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
            // reduce the degree (which is > 64) by dividing it by 64
            let factor = (-1i128).pow((deg / 64) as u32); // TODO too big ints?
            let new_deg = deg % 64;

            let mut term_poly_data = vec![Z_q::zero(); deg+1];
            term_poly_data[new_deg] = term*factor;
            let term_poly : R_q = R_q::new(term_poly_data);
            prod_rq = &prod_rq + &term_poly;
        }
        prod_rq
    }
}

impl std::ops::Mul<R_q> for &R_q {

    type Output = R_q;
        
    fn mul(self, rhs: R_q) -> Self::Output {
        // Custom multiplication logic.
        // We need to reduce by (X^d+1)

        // TODO add NTT logic here later.
        let prod : Polynomial<Z_q> = &self.0 * &rhs.0;
        let data_len = (&prod).data().len()+1;

        let mut prod_rq = R_q(prod);

        // REDUCE all terms of deg > 64
        for deg in (D as usize)+1..data_len {
            let term : Z_q = (&prod_rq).0.data()[deg].clone();
            // reduce the degree (which is > 64) by dividing it by 64
            let factor = (-1i128).pow((deg / 64) as u32); // TODO too big ints?
            let new_deg = deg % 64;

            let mut term_poly_data = vec![Z_q::zero(); deg+1];
            term_poly_data[new_deg] = term*factor;
            let term_poly : R_q = R_q::new(term_poly_data);
            prod_rq = &prod_rq + &term_poly;
        }
        prod_rq
    }
}

impl PartialEq for R_q {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}




impl fmt::Display for R_q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0.pretty("x"))
    }
}
