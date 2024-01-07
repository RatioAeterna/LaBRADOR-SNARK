use lazy_static::lazy_static;


// modulus of the ring of integers
pub const Q: i64 = 1 << 32; // 2^32, as used in the paper

// polynomial degree modulus
pub const D: i64 = 64; // used in the paper


// commitment ranks... TODO I have NO idea what these should actually be sized at. It doesn't seem
// to be described in the paper.
pub const KAPPA : i64 = 1024;
pub const KAPPA_1 : i64 = 1024;
pub const KAPPA_2 : i64 = 1024;



// number of functions 'f' in the family F of the principal relation R
pub const K : usize = 1;
// number of functions 'f-prime' in the family F' of the principal relation R
// TODO this MOST LIKELY SHOULD NOT BE CONSTANT, since it's based on the cardinality of F', which
// is whatever that happens to be.
// TODO actually, on second thought... this seems constant since it's not "whatever it happens to
// be" but whatever you set it to be.
pub const L : usize = 1;


// matrix dimensions: totally random for now and should be changed later
pub const N : usize = 128; // NUMBER of ROWS of S
pub const R : usize = 20; // Number of COLUMNS of S

// setting bound for SIS 2-norm
pub const BETA_BOUND : i64 = 65536;



// used for "constants" which cannot be evaluated at compile time
// most of these are the decomposition parameters described in 5.4
lazy_static! {
    // standard deviuation of the Z_q coefficients of the vectors s_i..
    // referred to as Gothic script s in section 5.4
    static ref STD : f64 = BETA_BOUND / ((R*N*D).sqrt());

    static ref B : i64 = ((12*R*TAU.sqrt()*STD).sqrt()).round();

    static ref T_1 : i64 = (Q.log10() / B.log10()).round();

    static ref B_1 : i64 = Q.powf((1.0 / T_1) as f64);

    static ref T_2 : i64 = (((24*N*D).sqrt()*(STD.pow(2))).log10() / B.log10()).round(); 

    static ref B_2 : i64 = ( ((25*N*D).sqrt()*(STD.pow(2))).powf((1.0 / T_2) as f64) ).round();

    static ref GAMMA : f64 = BETA_BOUND * TAU.sqrt();

    static ref GAMMA_1 : f64 = (((B_1.pow(2)*T_1) / 12.0)*R*KAPPA*D  + ((B_1.pow(2)*T_1) / 12.0)*((R.pow(2)+R)/2.0)*D).sqrt();

    static ref GAMMA_2 : f64 = (((B_1.pow(2)*T_1) / 12.0)*((R.pow(2)+R)/2.0)*D).sqrt();

    static ref BETA_PRIME : f64 = (((2.0/B.pow(2)) as f64)*GAMMA.pow(2) + GAMMA_1.pow(2) + GAMMA_2.pow(2)).sqrt();
}








