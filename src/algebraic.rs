use std::ops::Mul;
use std::fmt;
use ndarray::{Array2, Ix2, concatenate, Axis}; 
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use crate::constants::*;
use nalgebra::min;
use num_traits::Zero;
use num_traits::One;
use std::cmp::Ordering;
//use concrete_ntt::native64::Plan;

/*
 * Crate for various algebraic structure implementations
 * and related functions
 *
 */

pub fn mod_positive(dividend: i128, divisor: i128) -> i128 {
    let remainder = dividend % divisor;
    if remainder < 0 {
        remainder + divisor
    } else {
        remainder
    }
}


// The ring of integers modulo q (Q in constants.rs)
#[derive(Clone, Copy, Debug, Eq, Ord, PartialOrd)]
pub struct Z_q {
    value: i128,
}

impl Z_q {
    pub fn new(value: i128) -> Self {
        Z_q { value : mod_positive(value, Q) }
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
        self.value = mod_positive((self.value + other.value), Q);
    }
}

impl std::ops::SubAssign<Z_q> for Z_q {
    fn sub_assign(&mut self, other: Z_q) {
        self.value = mod_positive((self.value - other.value), Q);
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

impl From<Z_q> for i128 {
    fn from(item: Z_q) -> Self {
        item.value as i128
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
        self.value == mod_positive(*other, Q)
    }
}

impl PartialOrd<i128> for Z_q {
    fn partial_cmp(&self, other: &i128) -> Option<Ordering> {
        self.value.partial_cmp(&(mod_positive(*other, Q)))
    }
}



// Wrapper for Polynomial<Z_q> which specifically uses
// properties of R_q, aka Z[q]/(X^d+1)
#[derive(Clone, Debug)]
pub struct R_q(Polynomial<Z_q>);

impl R_q {
    // TODO maybe overload so you can pass in refs instead..
    pub fn new(coefficients: Vec<Z_q>) -> Self {
        R_q::reduction(coefficients)
    }

    pub fn data_vec(&self) -> Vec<Z_q> {
        self.0.data().to_vec()
    }

    pub fn eval(&self,  x: Z_q) -> Z_q {
        self.0.eval(x)
    }

    pub fn get_term_of_deg(&self, deg : usize) -> R_q {
        let datavec = self.0.data().to_vec();
        if (datavec.len() <= deg) { return R_q::zero(); }
        let coeff = datavec[deg];
        let mut ret_rq = vec![Z_q::zero(); deg+1];
        ret_rq[deg] = coeff;
        R_q::new(ret_rq)
    }

    fn reduction(coefficients: Vec<Z_q>) -> R_q {
        // REDUCE all terms of deg >= 64
        // form a new polynomial from the "valid" slice of data that we have, i.e., indices 0..D
        let data_len = coefficients.len();
        let valid_coeffs = coefficients[..min(data_len, 64)].to_vec();
        let sliced_poly : Polynomial<Z_q> = Polynomial::new(valid_coeffs);
        let mut reduced_rq : R_q = R_q(sliced_poly);
        if (data_len <= 64) {
            return reduced_rq;
        }

        for deg in (D as usize)..data_len {
            let term : Z_q = (&coefficients)[deg].clone();
            // reduce the degree (which is > 64) by dividing it by 64
            let factor = (-1i128).pow((deg / 64) as u32); // TODO too big ints?
            let new_deg = deg % 64;

            let mut term_poly_data = vec![Z_q::zero(); new_deg+1];
            term_poly_data[new_deg] = term*factor;
            let term_poly : R_q = R_q::new(term_poly_data);
            reduced_rq = &reduced_rq + &term_poly;
        }
        reduced_rq
    }

    // base multiplication function called from everywhere else
    fn multiply(lhs: &R_q, rhs: &R_q) -> R_q {
        // Custom multiplication logic.
        // We need to reduce by (X^d+1)
        // TODO add NTT logic here later.
        
        /*
        if(*NTT_ENABLED) {
            const N: usize = 32;
            let plan = Plan::try_new(D as usize).unwrap();

            let lhs_data = transform_slice_zq_to_u64((&lhs.0).data());
            let rhs_data = transform_slice_zq_to_u64((&rhs.0).data());
            let prod : Vec<u64> = vec![];

            
            negacylic_polymul(&prod, &lhs_data, &rhs_data);
        }
        */



        let prod : Polynomial<Z_q> = &lhs.0 * &rhs.0;
        R_q::new(prod.data().to_vec())
    }
}


/*
fn transform_slice_zq_to_u64(slice: &[Z_q]) -> Vec<u64> {
    slice.iter().map(|z| z.to_u128()).collect()
}
*/





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
        R_q::multiply(&self, &rhs) 
    }
}

impl std::ops::Mul<&R_q> for R_q {
    type Output = Self;

    fn mul(self, rhs: &R_q) -> Self::Output {
        R_q::multiply(&self, rhs) 
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
        R_q::multiply(self, &rhs)
    }
}

impl std::ops::Mul<R_q> for &R_q {

    type Output = R_q;
        
    fn mul(self, rhs: R_q) -> Self::Output {
        R_q::multiply(self, &rhs)
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
