use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;

pub struct Transcript {
    // fields (see protocol)
    z : Vec<Polynomial<i64>>,
}


pub struct State {
    // fields (see protocol)
}
