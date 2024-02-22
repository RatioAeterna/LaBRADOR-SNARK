use lazy_static::lazy_static;


// modulus of the ring of integers
pub const Q: i128 = 1 << 32; // 2^32, as used in the paper

// polynomial degree modulus
pub const D: i128 = 64; // used in the paper


// commitment ranks... TODO I have NO idea what these should actually be sized at. It doesn't seem
// to be described in the paper. Tried 1024, seems to be... a bit large. We'll see.
pub const KAPPA : usize = 128;
pub const KAPPA_1 : usize = 128;
pub const KAPPA_2 : usize = 128;

// used for challenge polynomial generation
pub const TAU : f64 = 71.0;
pub const T : f64 = 15.0;


// number of functions 'f' in the family F of the principal relation R
pub const K : usize = 1;
// number of functions 'f-prime' in the family F' of the principal relation R
// TODO this MOST LIKELY SHOULD NOT BE CONSTANT, since it's based on the cardinality of F', which
// is whatever that happens to be.
// TODO actually, on second thought... this seems constant since it's not "whatever it happens to
// be" but whatever you set it to be.
pub const L : usize = 1;


// matrix dimensions: totally random for now and should be changed later
//pub const N : usize = 128; // NUMBER of ROWS of S
pub const N : usize = 32; // NUMBER of ROWS of S
//pub const R : usize = 20; // Number of COLUMNS of S
pub const R : usize = 32; // Number of COLUMNS of S

// setting bound for SIS 2-norm
pub const BETA_BOUND : i128 = 65536;


// TODO it might be easier to read if you broke all of these down into functions e.g., 
// B_1: i128 = calculate_b1(Q, &T_1);

// used for "constants" which cannot be evaluated at compile time
// most of these are the decomposition parameters described in 5.4
lazy_static! {
    // standard deviuation of the Z_q coefficients of the vectors s_i..
    // referred to as Gothic script s in section 5.4
    pub static ref STD : f64 = (BETA_BOUND as f64) / ((((R*N) as f64)*D as f64).sqrt());

    pub static ref B : i128 = (((((12.*(R as f64)*TAU).sqrt()) as f64)*(*STD)).sqrt()).round() as i128;

    pub static ref T_1 : i128 = ((Q as f64).log10() / (*B as f64).log10()).round() as i128;

    // TODO do we need to round this instead? Unsure, will just truncate for now
    pub static ref B_1 : i128 = (Q as f64).powf((1.0 / (*T_1 as f64)) as f64) as i128;

    pub static ref T_2 : i128 = (((24.*((N as i128)*D) as f64).sqrt()*((*STD).powi(2))).log10() / (*B as f64).log10()).round() as i128; 

    pub static ref B_2 : i128 = ((((25*((N as i128)*D)) as f64).sqrt()*(STD.powi(2))).powf(1.0 / (*T_2 as f64))).round() as i128;

    pub static ref GAMMA : f64 = (BETA_BOUND as f64) * TAU.sqrt();

    pub static ref GAMMA_1 : f64 = ((((*B_1 as f64).powi(2)*(*T_1 as f64)) / 12.0)*(R as f64)*(KAPPA as f64)*(D as f64) + (((*B_1 as f64).powi(2)*(*T_1 as f64)) / 12.0)*((((R as f64).powi(2)+(R as f64)) as f64)/2.0)*(D as f64)).sqrt();

    pub static ref GAMMA_2 : f64 = ((((*B_1 as f64).powi(2)*(*T_1 as f64)) / 12.0)*(((R as f64).powi(2)+(R as f64))/2.0)*(D as f64)).sqrt();

    pub static ref BETA_PRIME : f64 = (((2.0/(*B as f64).powi(2)) as f64)*(*GAMMA as f64).powi(2) + (*GAMMA_1 as f64).powi(2) + (*GAMMA_2 as f64).powi(2)).sqrt();

    pub static ref NTT_ENABLED: bool = true;

    pub static ref AVX512_ENABLED: bool = false;
}








