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


