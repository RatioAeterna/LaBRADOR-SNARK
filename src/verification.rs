use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;

/*
pub fn valid_projection() : bool {


    compute_norm(&DVector::from_vec(projected_vec)) <= sqrt(128)*BETA_BOUND
}
*/
