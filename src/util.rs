use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;
// TODO yes, we're using this for now. No, it's not great. Yes, we will remove this later.
use std::clone::Clone;

use crate::constants::*;

// generates an element of Z_q[X]/(X^d+1)
// TODO: more sophisticated / cleaner random polynomial sampling
pub fn generate_polynomial(q : i64, d : i64) -> Polynomial<i64> {
    let max_degree = d-1; // TODO is this REALLY the max degree? I think NOT. Fix later.
    let mut rng = rand::thread_rng();
    //let num_terms = rng.gen_range(1, max_degree + 2);
    let mut coefficient_vec: Vec<i64> = Vec::new();
    for degree in 0..max_degree {
        coefficient_vec.push(rng.gen_range(0..q));
    }
    let mut poly = Polynomial::new(coefficient_vec);
    poly
}


pub fn generate_random_matrix(rows : i64, cols : i64, q : i64, d : i64) -> Array2<Polynomial<i64>> {
    let mut mat = Array2::from_elem((rows as usize,cols as usize), Polynomial::new(vec![])); 
    for i in 0..rows {
        for j in 0..cols {
            let poly = generate_polynomial(q,d);
            mat[[i as usize,j as usize]] = poly.clone();
        }
    }
    mat
}




// TODO yes we currently have to clone coeff_dist in when we want to re-use it, but it's NBD
// because it's usually very small. We don't want to use a mutable ref because that would be weird,
// since we don't want to change the underlying distribution data / remove elements if we're gonna re-use it later...
pub fn generate_polynomial_picky(q : i64, d : usize, mut coeff_dist : Vec<i64>) -> Polynomial<i64> {
    assert!(coeff_dist.len() == d, "Must have one coefficient for each degree of polynomial");
    let mut rng = rand::thread_rng();
    let mut coefficient_vec: Vec<i64> = Vec::new();
    for degree in 0..d {
        let random_index : usize = *(coeff_dist.as_slice().choose(&mut rng).unwrap()) as usize;
        let coeff = coeff_dist[random_index];
        coeff_dist.remove(random_index);
        let signed : bool = rng.gen();
        if (coeff > 0) && signed {
            coefficient_vec.push((-1*coeff) as i64);
        }
        else {
            coefficient_vec.push(coeff);
        }
    }
    let mut poly = Polynomial::new(coefficient_vec);
    poly
}

// The "Conjugation automorphism"
// This basically replaces all powers of X^n with X^{-n} and then reduces modulo X^d+1 (R's
// modulus) 
pub fn sigma_inv(a : &Polynomial<i64>) -> Polynomial<i64> {

    let poly_data_vec : Vec<i64> = a.data().to_vec();

    let mut new_coeff_vec : Vec<i64> = vec![0; 64]; 

    // terms of degree 0 unchanged
    // degree 1-64, we simply transform X^n to -X^{64-n}.

    for deg in 0..poly_data_vec.len() {
        let coeff = poly_data_vec[deg];
        if deg == 0 {
            new_coeff_vec[0] = coeff;
        }
        else {
            let new_deg = 64 - deg;
            new_coeff_vec[new_deg] = -1*coeff;
        }
    }
    let res_polynomial = Polynomial::new(new_coeff_vec);
    res_polynomial
}



pub fn multiply_poly_vec_ints(p_vec : &Vec<Polynomial<i64>>, ints: &Vec<i64>) -> Vec<Polynomial<i64>> {
    let mut p_vec_res : Vec<Polynomial<i64>> = vec![];

    for p in p_vec {
        p_vec_res.push(multiply_poly_ints(p, ints));
    }
    p_vec_res
}

pub fn multiply_poly_ints(p : &Polynomial<i64>, ints: &Vec<i64>) -> Polynomial<i64> {
    let mut p_res : Polynomial<i64> = Polynomial::new(vec![]);

    for coeff in ints {
        p_res = p_res + scale_polynomial(p, *coeff as f32);
    }
    p_res
}


// randomly samples n integers mod q, returns them as a vec
pub fn random_sample_Z_q(n : i64, q: i64) -> Vec<i64> {

    let mut rng = thread_rng();
    let dist = Uniform::from(0..q);

    let sample: Vec<i64> = (0..n).map(|_| rng.sample(dist)).collect();
    sample
}


// scales a given polynomial by a scale factor (usually < 1)
pub fn scale_polynomial(p : &Polynomial<i64>, s : f32) -> Polynomial<i64> {
    //let constant_polynomial : Polynomial<f32> = Polynomial::new(vec![s]);
    let poly_vec : Vec<i64> = p.data().to_vec(); 
    let scaled_poly_vec : Vec<i64> = poly_vec.iter().map(|&x| (((x as f32) * s).floor() as i64)).collect();
    return Polynomial::new(scaled_poly_vec);
}

// takes the 2-norm of a given polynomial (squared)
pub fn poly_norm(p : &Polynomial<i64>) -> f64 {
    let norm : i64 = p.data().to_vec().iter().map(|&x| x*x).sum::<i64>();
    norm as f64
}

// Does NOT square the norm
pub fn poly_norm_strict(p : &Polynomial<i64>) -> f64 {
    let strict_norm : f64 = ((p.data().to_vec().iter().map(|&x| x*x).sum::<i64>()) as f64).sqrt();
    strict_norm
}


// Does NOT square the norm. Integer vec version
pub fn l2_norm(vec : &Vec<i64>) -> f64 {
    let sum_of_squares: i64 = vec.iter().map(|&x| x * x).sum();
    (sum_of_squares as f64).sqrt()
}

// TODO we're using a statistical estimate of the sup to accomplish this, not sure if that's valid.
// should be since there's basically no other way to feasibly do this (I think). Check this.
pub fn operator_norm(c : &Polynomial<i64>) -> f64 {
    // here's the play.
    // compute n_samples sample polynomials r (where n_samples is some large number, but not too large)
    // compute the described ratio / inner product
    // keep a running estimate of the supremum. For sufficiently large N, it will be pretty approximate, and we can return this.
    // TODO tweak this value as needed
    let n_samples = 1000;

    let mut sup_estimate : f64 = 0.;

    for n in 0..n_samples {
        let r = generate_polynomial(Q, D);
        let norm : f64 = poly_norm_strict(&(c * &r));        
        let ratio : f64 = (norm / poly_norm_strict(&r)); 
        if ratio > sup_estimate {
            sup_estimate = ratio;
        }
    }
    sup_estimate
}

pub fn compute_total_norm(projection: &Array2<Polynomial<i64>>) -> f64 {
    let mut total_norm_squared: f64 = 0.0;
    for row in projection.outer_iter() {
        for poly in row.iter() {
            let norm = poly_norm(poly);
            total_norm_squared += norm;
        }
    }
    return f64::sqrt(total_norm_squared);
}

pub fn vec_inner_product(v1: &Vec<i64>, v2: &Vec<i64>) -> i64 {
    assert!(v1.len() == v2.len(), "inner product not defined on vectors of unequal length");
    let mut result : i64 = 0;
    for i in 0..v1.len() {
        let product = &v1[i] * &v2[i];
        result += product;
    }
    result
}


// elementwise multiplication of a vector of polynomials by a single polynomial. NOT an inner
// product.
pub fn poly_by_poly_vec(poly: &Polynomial<i64>, v1: &Vec<Polynomial<i64>>) -> Vec<Polynomial<i64>> {
    let mut v2 : Vec<Polynomial<i64>> = vec![];

    for i in 0..v1.len() {
        v2.push(poly * &v1[i]);
    }
    v2
}


pub fn add_poly_vec(v1: &Vec<Polynomial<i64>>, v2: &Vec<Polynomial<i64>>) -> Vec<Polynomial<i64>> {
    assert!(v1.len() == v2.len(), "summation not defined on vectors of unequal length");
    let mut v_res : Vec<Polynomial<i64>> = vec![];
    for i in 0..v1.len() {
        let new_p = &v1[i] + &v2[i];
        v_res.push(new_p);
    }
    v_res
}

pub fn add_vecs(v1: &Vec<i64>, v2: &Vec<i64>) -> Vec<i64> {
    assert!(v1.len() == v2.len(), "summation not defined on vectors of unequal length");
    let mut v_res : Vec<i64> = vec![];
    for i in 0..v1.len() {
        let new_p = &v1[i] + &v2[i];
        v_res.push(new_p);
    }
    v_res
}



/* Rarely, we need to do elementwise addition of a single polynomial to an entire vector of
 * polynomials.. kind of strange, but how this is described in the protocol */
pub fn add_poly_vec_by_poly(v1: &Vec<Polynomial<i64>>, p: &Polynomial<i64>) -> Vec<Polynomial<i64>> {
    let mut v_res : Vec<Polynomial<i64>> = vec![];
    for i in 0..v1.len() {
        let new_p = &v1[i] + p;
        v_res.push(new_p);
    }
    v_res
}





pub fn gen_empty_poly_vec(n : usize) -> Vec<Polynomial<i64>> {
    let mut v_res : Vec<Polynomial<i64>> = vec![];
    for i in 0..n {
        v_res.push(Polynomial::new(vec![]));
    }

    v_res
}

// used specifically do decompose e.g., the entirety of t_i = t_i^(0) + ... + t_i^(t_1-1)b_1^(t_1-1)
pub fn decompose_polynomial_vec(vec : &Vec<Polynomial<i64>>, base : i64, exp: i64) -> Vec<Vec<Polynomial<i64>>> {

    // TODO KAPPA is hardcoded here for now. Fix later.
    let mut res = vec![vec![Polynomial::new(vec![]); KAPPA as usize]; exp as usize];

    for i in 0..vec.len() {
        let t_i = &vec[i]; // polynomial we want to decompose
        let dec_t_i = decompose_polynomial(t_i, base, exp);
        for j in 0..(exp as usize) {
            // TODO yes, we're cloning here for now... we'll find a way to fix later.
            // Not too egregious in this instance. 
            res[j][i] = dec_t_i[j].clone(); 
        }
    }
    res
}

// computes the centered representative of a polynomial coefficient wrt a given base, i.e.,
// returns a value in the range of [-b/2, b/2]
pub fn centered_rep(mut val : i64, b: i64) -> i64 {
    if val > b / 2 {
        val -= b;
    }
    else if val <= (-b / 2) {
        val += b;
    }
    val
}


pub fn decompose_polynomial(p : &Polynomial<i64>, base : i64, exp: i64) -> Vec<Polynomial<i64>> {

    let poly_data_vec : Vec<i64> = p.data().to_vec();

    // this takes the form of, e.g., g_{ij} = g_{ij}^(0) + ... + g_{ij}^{t_2-1}b_2^{t_2-1}
    let mut decomp_vec : Vec<Polynomial<i64>> = vec![Polynomial::<i64>::new(vec![]); exp as usize];

    // we decompose each coefficient a_j of the polynomial (which we can express as a base K 

    for deg in 0..poly_data_vec.len() {
        let mut a_j = poly_data_vec[deg];
        while a_j != 0 {
            let a_jk = centered_rep(a_j % base, base);

            let mut a_jk_poly_data = vec![0; deg+1];
            a_jk_poly_data[deg] = a_jk;
            let a_jk_poly = Polynomial::new(a_jk_poly_data);
            decomp_vec[deg] = &decomp_vec[deg] + &a_jk_poly;
            a_j -= a_jk;
            a_j = (a_j / base); // this is the next coefficient to decompose
        }
    }
    decomp_vec
}

// Used to convert a vector of polynomials \vec{\b{s}}_i to another vector
// which is a concatenation of all the coefficients 
pub fn witness_coeff_concat(vec: &Vec<Polynomial<i64>>) -> Vec<i64> {
    // NOTE: we assume that the polynomials in vec will be of degree D... 
    // But we want to return a vec of degree N*D (where N is the assumed length of vec)..
    // So we recognize that for a polynomial ring mod X^D+1, X^d congruent to -1... so we can
    // rewrite as such.
    let mut coeffs : Vec<i64> = vec![];
    for j in 0..vec.len() {
        let poly : &Polynomial<i64> = &vec[j];
        let poly_vec_data : Vec<i64> = poly.data().to_vec();
        
        // TODO can we simplify this to poly_vec_data.len()? I don't think so, since data might
        // store fewer coeffs
        for deg in 0..(D as usize) {
            if deg < poly_vec_data.len() {
                // The case where we need to convert X^d -> -1 and add that to the constant term
                if (deg == 0) && (poly_vec_data.len() == (D as usize)+1) {
                    let constant_term = poly_vec_data[0]; 
                    let highest_term = poly_vec_data[(D as usize)];
                    let mut new_term = highest_term + constant_term;
                    // if negative, add Q until number is in [0, Q-1]
                    if(new_term < 0) {
                        while(new_term < 0) {
                            new_term += Q;
                        }
                    }
                    coeffs.push(new_term);
                }
                else {
                    coeffs.push(poly_vec_data[deg]);
                }
            }
            else {
                coeffs.push(0);
            }
        }
    }
    coeffs
}

// Convert the one-dimensional Array to a two-dimensional Array with one column
pub fn vec_to_column_array<T>(vec: &Vec<T>) -> Array2<T>
where
    T: Clone,
{
    // Convert the Vec<T> to a one-dimensional "column vector" Array
    // args to this function: num rows, num columns, vec to transform.
    let a = Array2::from_shape_vec((vec.len(), 1), vec.to_vec()).unwrap();
    a
}

// compute the dot product between two lists of polynomials
pub fn polynomial_vec_inner_product(v1: &[Polynomial<i64>], v2: &[Polynomial<i64>]) -> Polynomial<i64> {
    assert!(v1.len() == v2.len(), "inner product not defined on vectors of unequal length");
    let mut result = Polynomial::new(vec![]);
    for i in 0..v1.len() {
        let product = &v1[i] * &v2[i];
        result = result + product;
    }
    result
}

pub fn matmul(m1: &Array2<i64>, m2: &Array2<i64>) -> Array2<i64> {
    assert!(m1.shape()[1] == m2.shape()[0], "matrix product not defined on matrices which do not have dimension (m,n) x (n,k), for some m,n,k");
    let m = m1.shape()[0];
    let k = m2.shape()[1];

    let mut result = Array2::from_elem(Ix2(m, k), 0); 
    for i in 0..m {
        for j in 0..k {
            let v1 = m1.row(i).to_vec();
            let v2 = m1.column(j).to_vec();
            result[[i,j]] = vec_inner_product(&v1, &v2);
        }
    }
    result
}

// compute the matrix product between two polynomial MATRICES (Array2, slightly different)
pub fn polynomial_matrix_product(m1: &Array2<Polynomial<i64>>, m2: &Array2<Polynomial<i64>>) -> Array2<Polynomial<i64>> {
    assert!(m1.shape()[1] == m2.shape()[0], "matrix product not defined on matrices which do not have dimension (m,n) x (n,k), for some m,n,k");
    let m = m1.shape()[0];
    let k = m2.shape()[1];

    let mut result = Array2::from_elem(Ix2(m, k), Polynomial::new(vec![])); 
    for i in 0..m {
        for j in 0..k {
            let v1 = m1.row(i).to_vec();
            let v2 = m1.column(j).to_vec();
            result[[i,j]] = polynomial_vec_inner_product(&v1, &v2);
        }
    }
    result
}
