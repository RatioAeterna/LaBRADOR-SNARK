// modulus of the ring of integers
pub const Q: i64 = 128;

// polynomial degree modulus
pub const D: i64 = 5;


// number of functions 'f' in the family F of the principal relation R
pub const K : usize = 1;
// number of functions 'f-prime' in the family F' of the principal relation R
// TODO this MOST LIKELY SHOULD NOT BE CONSTANT, since it's based on the cardinality of F', which
// is whatever that happens to be.
pub const L : usize = 1;


// matrix dimensions: totally random for now and should be changed later
pub const N : usize = 128; // NUMBER of ROWS of S
pub const R : usize = 20; // Number of COLUMNS of S

// setting bound for SIS 2-norm
pub const BETA_BOUND : i64 = 65536;
