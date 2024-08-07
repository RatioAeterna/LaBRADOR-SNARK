use crate::util::*;
use concrete_ntt::native64::Plan32;
use lazy_static::lazy_static;
use num_prime::nt_funcs::is_prime;
use once_cell::sync::OnceCell;
use std::sync::atomic::AtomicBool;

pub static VERBOSE: OnceCell<bool> = OnceCell::new();

pub fn is_verbose() -> bool {
    *VERBOSE.get().unwrap_or(&false)
}

// polynomial degree modulus
pub const D: usize = 64; // used in the paper
//pub const D: usize = 2; 


// used for challenge polynomial generation
// TODO currently fixed values to work with D=64, Q around 2^32, but we likely need
// different values for different ring dimension, etc.
pub const TAU: f64 = 71.0;
pub const T: f64 = 15.0;

// number of functions 'f' in the family F of the principal relation R
pub const K: usize = 1;
// number of functions 'f-prime' in the family F' of the principal relation R
// TODO this MOST LIKELY SHOULD NOT BE CONSTANT, since it's based on the cardinality of F', which
// is whatever that happens to be.
// TODO actually, on second thought... this seems constahe cost of precision. An f64 can represent numbers much larger than an i128, but it does so with a degree of approximation for extremely large or small numbent since it's not "whatever it happens to
// be" but whatever you set it to be.
pub const L: usize = 1;


// NOTES here on what exactly constitutes the set "R_q" / how we properly choose Q
// for Z[x] / (x^64+1), where D=64...
//
// - take a generator of (Z / 128Z)^*
// - square it -> gives you an element of order 32 in (Z/128Z)^*
//
// want that g^2 + 128k to be prime in Z.. any prime of that form works

// TODO this generator finding code appears to be broken. For now, we'll just select the closest
// prime.
pub fn find_suitable_prime(start: i128) -> i128 {
    //println!("in here!");

    // first, find a generator of (Z / 2D Z)^*...

    /*
    let g = find_generator(2*(D as i128));
    println!("test 3");
    let base = g*g;
    
    println!("generator is... {}", g);

    // now, add 2Dk for k in Z until we get a prime... and total value is >= start
    let mut k = ((start - base) as f64 / (2.0 * D as f64)).ceil() as i128; // starting value
    
    // TODO slight bug where this is always guaranteed to be larger (fine) but k sometimes one
    // value higher than it needs to be to start.. rounding error
    let mut q = base + 2*(D as i128)*k;
    assert!(q >= start, "q not geq requested start value.. inequality bug");
    */
    let mut q = start;

    while q < i128::MAX {
        // TODO does having the 'is_suitable' constraint even help, if we know we're actually using
        // R_q?
        /*
        if is_prime(&(q as u128), None).probably() && is_suitable(q) {
            return q;
        }
        */
        if is_prime(&(q as u128), None).probably() {
            return q;
        }
        /*
        q += 2*(D as i128);
        k += 1;
        println!("k = {}, q = {}", k, q);
        */
        q += 1;
    }
    //println!("value for k: {}", k);
    panic!("Error: No suitable value for Q found!");
}

fn mod_exp(base: i128, exp: i128, modulus: i128) -> i128 {
    let mut result = 1;
    let mut base = base % modulus;
    let mut exp = exp;

    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        exp = exp >> 1;
        base = (base * base) % modulus;
    }
    result
}

fn is_generator(g: i128) -> bool {
    let phi = (D as i128); // Euler's totient function for 2D
    let mut factors = vec![];

    // Find factors of phi
    let mut i = 2;
    let mut n = phi;
    while i * i <= n {
        if n % i == 0 {
            factors.push(i);
            while n % i == 0 {
                n /= i;
            }
        }
        i += 1;
    }
    if n > 1 {
        factors.push(n);
    }

    println!("factors: {:?}", factors);

    // Check if g^(phi/factor) != 1 (mod 2D) for all prime factors of phi
    let modulus = 2 * (D as i128);
    for factor in factors {
        println!("mod_exp: {}", mod_exp(g, phi / factor, modulus));
        if mod_exp(g, phi / factor, modulus) == 1 {
            println!("returning false");
            return false;
        }
    }
    true
}




/*
fn is_generator(g: i128) -> bool {
    let phi = (D as i128); // Euler's totient function for 2D
    let mut factors = vec![];

    // Find factors of phi
    let mut i = 2;
    let mut n = phi;
    while i * i <= n {
        if n % i == 0 {
            factors.push(i);
            while n % i == 0 {
                n /= i;
            }
        }
        i += 1;
    }
    if n > 1 {
        factors.push(n);
    }

    // Check if g^(phi/factor) != 1 (mod 2D) for all prime factors of phi
    for factor in factors {
        if g.pow((phi / factor) as u32) % (2*(D as i128)) == 1 {
            println!("returning false");
            return false;
        }
    }
    true
}
*/


// TODO refactor so this is 'find_odd_generator' and is prioritizing the parity check to save
// computation
fn find_generator(max_val: i128) -> i128 {
    println!("in here 2");
    for g in 1..max_val {
        println!("g=.. {}", g);
        if is_generator(g) && (g % 2 == 1){
            return g;
        }
    }
    0
}


// TODO EVENTUALLY we'll want to make sure it actually splits R_q into two irreducible ideals, etc.
pub fn is_suitable(candidiate_q: i128) -> bool {
    return mod_positive(candidiate_q, 2 * (D as i128)) == 1;
}

lazy_static! {
    // modulus of the ring of integers
    pub static ref Q: i128 = find_suitable_prime((1<<13)-1);
    //pub static ref Q: i128 = find_suitable_prime((1<<32)-1);
    pub static ref PLAN: Plan32 = Plan32::try_new(D as usize).unwrap();
}

pub static NTT_ENABLED: AtomicBool = AtomicBool::new(false);
pub static MOD_SUSPENSION: AtomicBool = AtomicBool::new(false);

// used for "constants" which cannot be evaluated at compile time
// many of these are the decomposition parameters described in 5.4
pub struct RuntimeConstants {
    pub N: usize, 
    pub R: usize, 

    // setting bound for SIS 2-norm
    // See Theorem 5.1 for description of upper bound..
    pub BETA_BOUND: i128,
    // standard deviation of the Z_q coefficients of the vectors s_i..
    // referred to as Gothic script s in section 5.4
    pub STD: f64,
    pub B: i128,
    pub T_1: i128,
    // TODO do we need to round this instead? Unsure, will just truncate for now
    pub B_1: i128,
    pub T_2: i128,
    pub B_2: i128,
    pub GAMMA: f64,
    pub GAMMA_1: f64,
    pub GAMMA_2: f64,
    pub BETA_PRIME: f64,

    // commitment ranks...
    pub KAPPA: usize,
    pub KAPPA_1: usize,
    pub KAPPA_2: usize,
}

impl RuntimeConstants {

    pub fn new(N: usize, R: usize) -> Self {
        // NOTE: commitment rank must be AT LEAST as large as the lattice dimension in order for
        // commitments to be binding.
        let KAPPA : usize = N*D;
        let KAPPA_1 : usize = N*D;
        let KAPPA_2 : usize = N*D;

        let BETA_BOUND : i128 = ((30.0/128.0 as f64).sqrt()*(*Q as f64)/125.0).floor() as i128;
        let STD : f64 = (BETA_BOUND as f64) / ((R*N*D) as f64).sqrt();
        let B : i128 = (((((12.*(R as f64)*TAU).sqrt()) as f64)*(STD)).sqrt()).round() as i128;
        let T_1 : i128 = ((*Q as f64).log2() / (B as f64).log2()).round() as i128;
        let B_1 : i128 = (*Q as f64).powf((1.0 / (T_1 as f64)) as f64) as i128;
        let T_2 : i128 = (((24.*(N*D) as f64).sqrt()*((STD).powi(2))).log2() / (B as f64).log2()).round() as i128;
        let B_2 : i128 = ((((24*(N*D)) as f64).sqrt()*(STD.powi(2))).powf(1.0 / (T_2 as f64))).round() as i128;

        // NOTE.. all of these are SQUARED gamma, gamma_1, etc.
        let GAMMA : f64 = ((BETA_BOUND*BETA_BOUND) as f64) * TAU;
        let GAMMA_1 : f64 = (((B_1 as f64).powi(2)*(T_1 as f64)) / 12.0)*(R as f64)*(KAPPA as f64)*(D as f64) + (((B_2 as f64).powi(2)*(T_2 as f64)) / 12.0)*((((R as f64).powi(2)+(R as f64)) as f64)/2.0)*(D as f64);
        let GAMMA_2 : f64 = (((B_1 as f64).powi(2)*(T_1 as f64)) / 12.0)*(((R as f64).powi(2)+(R as f64))/2.0)*(D as f64);
        // NOTE... this is ACTUALLY beta prime squared!
        let BETA_PRIME : f64 = (((2.0/(B as f64).powi(2)) as f64)*GAMMA + GAMMA_1 + GAMMA_2);
        /*
        let GAMMA : f64 = (BETA_BOUND as f64) * TAU.sqrt();
        let GAMMA_1 : f64 = ((((B_1 as f64).powi(2)*(T_1 as f64)) / 12.0)*(R as f64)*(KAPPA as f64)*(D as f64) + (((B_2 as f64).powi(2)*(T_2 as f64)) / 12.0)*((((R as f64).powi(2)+(R as f64)) as f64)/2.0)*(D as f64)).sqrt();
        let GAMMA_2 : f64 = ((((B_1 as f64).powi(2)*(T_1 as f64)) / 12.0)*(((R as f64).powi(2)+(R as f64))/2.0)*(D as f64)).sqrt();
        let BETA_PRIME : f64 = (((2.0/(B as f64).powi(2)) as f64)*(GAMMA as f64).powi(2) + (GAMMA_1 as f64).powi(2) + (GAMMA_2 as f64).powi(2)).sqrt();

        */

        Self { N, R, BETA_BOUND, STD, B, T_1, B_1, T_2, B_2, GAMMA, GAMMA_1, GAMMA_2, BETA_PRIME, KAPPA, KAPPA_1, KAPPA_2 }
    }
}
