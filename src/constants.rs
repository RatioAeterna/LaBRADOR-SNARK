use crate::util::*;
use lazy_static::lazy_static;
use std::sync::atomic::{AtomicBool, Ordering};
use concrete_ntt::native64::Plan32;
use num_prime::{PrimalityTestConfig};
use num_prime::nt_funcs::{is_prime};



// polynomial degree modulus
//pub const D: i128 = 64; // used in the paper
pub const D: i128 = 1; // DEBUG


// commitment ranks... 
// TODO resize based on application of core-svp methodology corresponding to adequate security
// parameter
pub const KAPPA : usize = 2;
pub const KAPPA_1 : usize = 2;
pub const KAPPA_2 : usize = 2;
/*
pub const KAPPA : usize = 256;
pub const KAPPA_1 : usize = 256;
pub const KAPPA_2 : usize = 256;
*/

// used for challenge polynomial generation
pub const TAU : f64 = 71.0;
pub const T : f64 = 15.0;


// number of functions 'f' in the family F of the principal relation R
pub const K : usize = 1;
// number of functions 'f-prime' in the family F' of the principal relation R
// TODO this MOST LIKELY SHOULD NOT BE CONSTANT, since it's based on the cardinality of F', which
// is whatever that happens to be.
// TODO actually, on second thought... this seems constahe cost of precision. An f64 can represent numbers much larger than an i128, but it does so with a degree of approximation for extremely large or small numbent since it's not "whatever it happens to
// be" but whatever you set it to be.
pub const L : usize = 1;


// matrix dimensions: totally random for now and should be changed later
pub const N : usize = 2; // NUMBER of ROWS of S
//pub const N : usize = 2; // NUMBER of ROWS of S
pub const R : usize = 2; // Number of COLUMNS of S
//pub const R : usize = 2; // Number of COLUMNS of S



pub fn find_suitable_prime(start : i128) -> i128 {
    let mut q = start;

    while q < i128::MAX {
        if is_prime(&(q as u128), None).probably() && is_suitable(q) {
            return q;
        }
        q += 2;
    }
    panic!("Error: No suitable value for Q found!");
}

// TODO EVENTUALLY we'll want to make sure it actually splits R_q into two irreducible ideals, etc.
pub fn is_suitable(candidiate_q : i128) -> bool {
    return (mod_positive(candidiate_q, 2*D) == 1);
}


// TODO it might be easier to read if you broke all of these down into functions e.g., 
// B_1: i128 = calculate_b1(Q, &T_1);

// used for "constants" which cannot be evaluated at compile time
// most of these are the decomposition parameters described in 5.4
lazy_static! {

    // modulus of the ring of integers
    pub static ref Q: i128 = find_suitable_prime((1<<10)-1);
    //pub static ref Q: i128 = find_suitable_prime((1<<32)-1);

    // setting bound for SIS 2-norm
    // See Theorem 5.1 for description of upper bound..
    pub static ref BETA_BOUND : i128 = ((30.0/128.0 as f64).sqrt()*(*Q as f64)/125.0).floor() as i128;

    // standard deviuation of the Z_q coefficients of the vectors s_i..
    // referred to as Gothic script s in section 5.4
    pub static ref STD : f64 = (*BETA_BOUND as f64) / ((((R*N) as f64)*D as f64).sqrt());

    pub static ref B : i128 = (((((12.*(R as f64)*TAU).sqrt()) as f64)*(*STD)).sqrt()).round() as i128;

    pub static ref T_1 : i128 = ((*Q as f64).log2() / (*B as f64).log2()).round() as i128;

    // TODO do we need to round this instead? Unsure, will just truncate for now
    pub static ref B_1 : i128 = (*Q as f64).powf((1.0 / (*T_1 as f64)) as f64) as i128;

    pub static ref T_2 : i128 = (((24.*((N as i128)*D) as f64).sqrt()*((*STD).powi(2))).log2() / (*B as f64).log2()).round() as i128; 

    pub static ref B_2 : i128 = ((((25*((N as i128)*D)) as f64).sqrt()*(STD.powi(2))).powf(1.0 / (*T_2 as f64))).round() as i128;

    pub static ref GAMMA : f64 = (*BETA_BOUND as f64) * TAU.sqrt();

    pub static ref GAMMA_1 : f64 = ((((*B_1 as f64).powi(2)*(*T_1 as f64)) / 12.0)*(R as f64)*(KAPPA as f64)*(D as f64) + (((*B_1 as f64).powi(2)*(*T_1 as f64)) / 12.0)*((((R as f64).powi(2)+(R as f64)) as f64)/2.0)*(D as f64)).sqrt();

    pub static ref GAMMA_2 : f64 = ((((*B_1 as f64).powi(2)*(*T_1 as f64)) / 12.0)*(((R as f64).powi(2)+(R as f64))/2.0)*(D as f64)).sqrt();

    pub static ref BETA_PRIME : f64 = (((2.0/(*B as f64).powi(2)) as f64)*(*GAMMA as f64).powi(2) + (*GAMMA_1 as f64).powi(2) + (*GAMMA_2 as f64).powi(2)).sqrt();

    pub static ref PLAN : Plan32 = Plan32::try_new(D as usize).unwrap();


}

pub static NTT_ENABLED: AtomicBool = AtomicBool::new(false);

pub static MOD_SUSPENSION: AtomicBool = AtomicBool::new(false);
