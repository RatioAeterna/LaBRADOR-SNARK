use std::ops::Mul;
use std::fmt;
use ndarray::{Array2, Ix2, concatenate, Axis}; 
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use crate::constants::*;
use num_traits::Zero;

/*
 * Functionality for the proof's recursive step
 *
*/


