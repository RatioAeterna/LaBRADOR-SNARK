use std::ops::Mul;
use std::fmt;
use ndarray::{Array2, Ix2, concatenate, Axis}; 
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use crate::constants::*;
use num_traits::Zero;




/*
 * Various functionality for assessing performance of this construction,
 * Both in terms of runtime and proof size, and logging the results, etc.
*/

/*
pub fn recursive_step(proof: &Transcript) -> Transcript {

    


    // size of the new witness (reduced from R)
    let R_prime : usize = 2*upsilon + mu;

    let N_prime : usize = max(ceil(n/upsilon), ceil(m/mu));
    
    // size of the new target relation
    let K_prime : usize = KAPPA + KAPPA_1 + KAPPA_2 + 3;








}
*/







