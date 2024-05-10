use crate::constants::*;
use crate::util::*;
use nalgebra::min;
use num_traits::One;
use num_traits::Zero;
use polynomial::Polynomial;
use std::cmp::Ordering;
use std::fmt;
use std::sync::atomic::Ordering as AtomicOrdering;

//#[cfg(test)]
use proptest::prelude::*;

/*
 * Crate for various algebraic structure implementations
 * and related functions
 *
 */

// The ring of integers modulo q (Q in constants.rs)
#[derive(Clone, Copy, Debug, Eq, Ord, PartialOrd)]
pub struct Zq {
    value: i128,
}

impl Zq {
    pub fn new(value: i128) -> Self {
        if MOD_SUSPENSION.load(AtomicOrdering::SeqCst) {
            Zq { value }
        } else {
            Zq {
                value: mod_positive(value, *Q),
            }
        }
    }
    pub fn lift(vec: &Vec<i128>) -> Vec<Zq> {
        let mut res_vec: Vec<Zq> = vec![];

        for i in 0..vec.len() {
            res_vec.push(Zq::from(vec[i]));
        }
        res_vec
    }
    pub fn lift_inv(vec: &Vec<Zq>) -> Vec<i128> {
        let mut res_vec: Vec<i128> = vec![];
        for i in 0..vec.len() {
            res_vec.push(vec[i].value);
        }
        res_vec
    }
}

impl std::ops::Neg for Zq {
    type Output = Zq;
    fn neg(self) -> Zq {
        Zq {
            value: if self.value == 0 { 0 } else { *Q - self.value },
        }
    }
}

impl Zero for Zq {
    fn zero() -> Self {
        Zq { value: 0 }
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }
}

impl std::ops::Rem<i128> for Zq {
    type Output = Zq;

    fn rem(self, rhs: i128) -> Self::Output {
        Zq::new(self.value % rhs)
    }
}

impl One for Zq {
    fn one() -> Self {
        Zq { value: 1 }
    }
}

impl std::iter::Sum for Zq {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl std::ops::Add for Zq {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::new(self.value + other.value)
    }
}

impl std::ops::Add<i128> for Zq {
    type Output = Self;

    fn add(self, other: i128) -> Self::Output {
        Self::new(self.value + other)
    }
}

impl std::ops::Sub for Zq {
    type Output = Self;

    // TODO is this correct? Double check
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.value - other.value)
    }
}

impl std::ops::Mul for Zq {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        Self::new(self.value * other.value)
    }
}

impl std::ops::Mul<i128> for Zq {
    type Output = Self;

    fn mul(self, other: i128) -> Self::Output {
        Self::new(self.value * other)
    }
}

impl std::ops::Div<i128> for Zq {
    type Output = Self;

    fn div(self, rhs: i128) -> Self::Output {
        Self::new(self.value / rhs)
    }
}

impl std::ops::Div for Zq {
    type Output = Self;

    fn div(self, rhs: Zq) -> Self::Output {
        /*
        let rhs_val_bigint = BigInt::from(rhs.value);
        let q_bigint = BigInt::from(*Q);
        let rhs_inv : i128 = rhs_val_bigint.modpow(&(q_bigint.clone() - rhs_val_bigint.clone()), &q_bigint).to_i128().unwrap();
        let div_val = self.value * rhs_inv;
        Self::new(div_val)
        */
        Self::new(self.value / rhs.value)
    }
}

impl std::ops::Add for &Zq {
    type Output = Zq;

    fn add(self, other: Self) -> Self::Output {
        Zq::new(self.value + other.value)
    }
}

impl std::ops::Sub for &Zq {
    type Output = Zq;

    // TODO is this correct? Double check
    fn sub(self, other: Self) -> Self::Output {
        Zq::new(self.value - other.value)
    }
}

impl std::ops::Mul for &Zq {
    type Output = Zq;

    fn mul(self, other: Self) -> Self::Output {
        Zq::new(self.value * other.value)
    }
}

impl std::ops::AddAssign<Zq> for Zq {
    fn add_assign(&mut self, other: Zq) {
        self.value = mod_positive(self.value + other.value, *Q);
    }
}

impl std::ops::SubAssign<Zq> for Zq {
    fn sub_assign(&mut self, other: Zq) {
        self.value = mod_positive(self.value - other.value, *Q);
    }
}

impl fmt::Display for Zq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //println!("DISPLAYING (the z_q...) {}", self.value);
        write!(f, "{}", self.value)
    }
}

impl From<&Zq> for usize {
    fn from(item: &Zq) -> Self {
        item.value as usize
    }
}

impl From<&Zq> for u64 {
    fn from(item: &Zq) -> Self {
        item.value as u64
    }
}

impl From<Zq> for f32 {
    fn from(item: Zq) -> Self {
        item.value as f32
    }
}

impl From<&Zq> for f32 {
    fn from(item: &Zq) -> Self {
        item.value as f32
    }
}

impl From<Zq> for f64 {
    fn from(item: Zq) -> Self {
        item.value as f64
    }
}

impl From<Zq> for i128 {
    fn from(item: Zq) -> Self {
        item.value as i128
    }
}

impl From<f32> for Zq {
    fn from(item: f32) -> Self {
        Self::new(item as i128)
    }
}

impl From<f64> for Zq {
    fn from(item: f64) -> Self {
        Self::new(item as i128)
    }
}

impl From<i32> for Zq {
    fn from(item: i32) -> Self {
        Self::new(item as i128)
    }
}

impl From<i128> for Zq {
    fn from(item: i128) -> Self {
        Self::new(item)
    }
}

impl PartialEq for Zq {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl PartialEq<i128> for Zq {
    fn eq(&self, other: &i128) -> bool {
        self.value == mod_positive(*other, *Q)
    }
}

impl PartialOrd<i128> for Zq {
    fn partial_cmp(&self, other: &i128) -> Option<Ordering> {
        self.value.partial_cmp(&(mod_positive(*other, *Q)))
    }
}

// Wrapper for Polynomial<Zq> which specifically uses
// properties of Rq, aka Z[q]/(X^d+1)
#[derive(Clone, Debug)]
pub struct Rq(Polynomial<Zq>);

impl Rq {
    // TODO maybe overload so you can pass in refs instead..
    pub fn new(coefficients: Vec<Zq>) -> Self {
        if MOD_SUSPENSION.load(AtomicOrdering::SeqCst) {
            let poly: Polynomial<Zq> = Polynomial::new(coefficients);
            Rq(poly)
        } else {
            Rq::reduction(coefficients)
        }
    }

    // assumes that moduluses are re-enabled.
    pub fn recompute_mod(&self) -> Self {
        let mut datavec = self.data_vec();
        for i in 0..datavec.len() {
            datavec[i] = Zq::new(i128::from(datavec[i]));
        }
        Rq::reduction(datavec)
    }

    pub fn data_vec(&self) -> Vec<Zq> {
        self.0.data().to_vec()
    }

    // TODO using this is maybe slightly dangerous, since coefficients could potentially be larger
    // than max u64
    pub fn raw_coeffs(&self) -> Vec<u64> {
        let datavec = self.0.data().to_vec();
        datavec.into_iter().map(|x| x.value as u64).collect()
    }

    pub fn eval(&self, x: Zq) -> Zq {
        self.0.eval(x)
    }

    pub fn get_term_of_deg(&self, deg: usize) -> Rq {
        let datavec = self.0.data().to_vec();
        if datavec.len() <= deg {
            return Rq::zero();
        }
        let coeff = datavec[deg];
        let mut ret_rq = vec![Zq::zero(); deg + 1];
        ret_rq[deg] = coeff;
        Rq::new(ret_rq)
    }

    fn reduction(coefficients: Vec<Zq>) -> Rq {
        // REDUCE all terms of deg >= D
        // form a new polynomial from the "valid" slice of data that we have, i.e., indices 0..D
        let data_len = coefficients.len();
        let valid_coeffs = coefficients[..min(data_len, D as usize)].to_vec();
        // TODO fix this cloning
        let sliced_poly: Polynomial<Zq> = Polynomial::new(valid_coeffs.clone());
        let mut reduced_rq: Rq = Rq(sliced_poly);
        if data_len <= (D as usize) {
            //println!("valid coeffs... {:?}", valid_coeffs);
            //println!("reduced r_q... {}", reduced_rq);
            return reduced_rq;
        }

        for deg in (D as usize)..data_len {
            let term: Zq = (&coefficients)[deg].clone();
            // reduce the degree (which is > D) by dividing it by D
            let factor = (-1i128).pow((deg / D as usize) as u32); // TODO too big ints?
            let new_deg = deg % (D as usize);

            let mut term_poly_data = vec![Zq::zero(); new_deg + 1];
            term_poly_data[new_deg] = term * factor;
            let term_poly: Rq = Rq::new(term_poly_data);
            reduced_rq = &reduced_rq + &term_poly;
        }
        reduced_rq
    }

    // base multiplication function called from everywhere else
    fn multiply(lhs: &Rq, rhs: &Rq) -> Rq {
        // Custom multiplication logic.
        // We need to reduce by (X^d+1), which we do in Rq::new()
        // TODO add NTT logic here later.

        if NTT_ENABLED.load(AtomicOrdering::SeqCst) {
            let lhs_data = transform_slice_zq_to_u64((&lhs.0).data());
            let rhs_data = transform_slice_zq_to_u64((&rhs.0).data());
            let mut prod: Vec<u64> = vec![0; D as usize];

            (*PLAN).negacyclic_polymul(&mut prod, &lhs_data, &rhs_data);
        }

        let prod: Polynomial<Zq> = &lhs.0 * &rhs.0;
        Rq::new(prod.data().to_vec())
    }
}

fn transform_slice_zq_to_u64(slice: &[Zq]) -> Vec<u64> {
    slice.iter().map(|z| u64::from(z)).collect()
}

impl Zero for Rq {
    fn zero() -> Self {
        Rq::new(vec![])
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl std::ops::Add for Rq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let sum: Polynomial<Zq> = self.0 + rhs.0;
        let sum_rq = Rq(sum);
        sum_rq
    }
}

impl std::ops::Add<&Rq> for Rq {
    type Output = Self;

    fn add(self, rhs: &Rq) -> Self::Output {
        let sum: Polynomial<Zq> = self.0 + &rhs.0;
        let sum_rq = Rq(sum);
        sum_rq
    }
}

impl std::ops::Sub for Rq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let sum: Polynomial<Zq> = self.0 - rhs.0;
        let sum_rq = Rq(sum);
        sum_rq
    }
}

impl std::ops::Mul for Rq {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Rq::multiply(&self, &rhs)
    }
}

impl std::ops::Mul<&Rq> for Rq {
    type Output = Self;

    fn mul(self, rhs: &Rq) -> Self::Output {
        Rq::multiply(&self, rhs)
    }
}

impl std::ops::Add for &Rq {
    type Output = Rq;

    fn add(self, rhs: Self) -> Self::Output {
        let sum: Polynomial<Zq> = &self.0 + &rhs.0;
        let sum_rq = Rq(sum);
        sum_rq
    }
}

impl std::ops::Add<Rq> for &Rq {
    type Output = Rq;

    fn add(self, rhs: Rq) -> Self::Output {
        let sum: Polynomial<Zq> = &self.0 + rhs.0;
        let sum_rq = Rq(sum);
        sum_rq
    }
}

impl std::ops::Sub for &Rq {
    type Output = Rq;

    fn sub(self, rhs: Self) -> Self::Output {
        let minus: Polynomial<Zq> = &self.0 - &rhs.0;
        let minus_rq = Rq(minus);
        minus_rq
    }
}

impl std::ops::Mul for &Rq {
    type Output = Rq;

    fn mul(self, rhs: Self) -> Self::Output {
        Rq::multiply(self, &rhs)
    }
}

impl std::ops::Mul<Rq> for &Rq {
    type Output = Rq;

    fn mul(self, rhs: Rq) -> Self::Output {
        Rq::multiply(self, &rhs)
    }
}

impl PartialEq for Rq {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl fmt::Display for Rq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //println!("DISPLAYING!! :) {:?}", self.0.data().to_vec());
        write!(f, "{}", self.0.pretty("x"))
    }
}

// Proptest trait implementations...
//
//

//#[cfg(test)]
impl Arbitrary for Zq {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_args: Self::Parameters) -> Self::Strategy {
        any::<i128>().prop_map(|x| Zq::new(x)).boxed()
    }
}

//#[cfg(test)]
impl Arbitrary for Rq {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_args: Self::Parameters) -> Self::Strategy {
        prop::collection::vec(any::<Zq>(), D as usize)
            .prop_map(Rq::new)
            .boxed()
    }
}
