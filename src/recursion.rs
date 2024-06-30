use std::ops::Mul;
use std::fmt;
use ndarray::{Array2, Ix2, concatenate, Axis}; 
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use crate::constants::*;
use num_traits::Zero;


pub fn recursive_step(proof: &Transcript) -> Transcript {
    // size of the new witness (reduced from R)
    let r_prime : usize = 2*nu + mu;

    let n_prime : usize = max(ceil(n/nu), ceil(m/mu));
    
    // size of the new target relation
    let k_prime : usize = KAPPA + KAPPA_1 + KAPPA_2 + 3;

    // QUESTION: what happens to commitment ranks KAPPA on subsequent iterations?
    // - I think we recompute as KAPPA = n_prime*D

    // TODO how do we compute 'nu'?
    // TODO how do we compute 'mu'?
    // TODO how are we "pruning out" the rest of the witness, and of the target relation?







}







