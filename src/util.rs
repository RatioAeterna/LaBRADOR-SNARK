use ndarray::{Array2, Ix2, concatenate};
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::{Distribution, WeightedIndex, Uniform};
// TODO yes, we're using this for now. No, it's not great. Yes, we will remove this later.
use std::clone::Clone;
use num_traits::Zero;
use std::ops::{AddAssign, Mul};

use crate::algebraic::*;
use crate::constants::*;

// generates an elemen100 of Z_q[X]/(X^d+1)
// TODO: more sophisticated / cleaner random polynomial sampling
pub fn generate_polynomial(q : i128, d : i128) -> R_q {
    let mut rng = rand::thread_rng();
    //let num_terms = rng.gen_range(1, max_degree + 2);
    let mut coefficient_vec: Vec<Z_q> = Vec::new();
    for degree in 0..d {
        coefficient_vec.push(Z_q::from(rng.gen_range(0..q)));
    }
    let mut poly = R_q::new(coefficient_vec);
    poly
}

pub fn generate_sparse_polynomial(q : i128, d : i128) -> R_q {
    let mut rng = rand::thread_rng();
    //let num_terms = rng.gen_range(1, max_degree + 2);
    let mut coefficient_vec: Vec<Z_q> = Vec::new();
    for degree in 0..d {
        // 90% chance to be zero
        if rng.gen_bool(0.9) {
            coefficient_vec.push(Z_q::zero());
        }
        else {
            coefficient_vec.push(Z_q::from(rng.gen_range(0..q/10000000)));
        }
    }
    let mut poly = R_q::new(coefficient_vec);
    poly
}


pub fn generate_random_matrix(rows : usize, cols : usize, q : i128, d : i128) -> Array2<R_q> {
    let mut mat = Array2::from_elem((rows,cols), R_q::new(vec![])); 
    for i in 0..rows {
        for j in 0..cols {
            let poly = generate_polynomial(q,d);
            mat[[i,j]] = poly.clone();
        }
    }
    mat
}




// TODO yes we currently have to clone coeff_dist in when we want to re-use it, but it's NBD
// because it's usually very small. We don't want to use a mutable ref because that would be weird,
// since we don't want to change the underlying distribution data / remove elements if we're gonna re-use it later...
pub fn generate_polynomial_picky(q : i128, d : usize, mut coeff_dist : Vec<Z_q>) -> R_q {
    assert!(coeff_dist.len() == d, "Must have one coefficient for each degree of polynomial");
    let mut rng = rand::thread_rng();
    let mut coefficient_vec: Vec<Z_q> = Vec::new();
    for degree in 0..d {
        let random_index = rng.gen_range(0..coeff_dist.len());
        //let random_index : usize = usize::from(coeff_dist.as_slice().choose(&mut rng).unwrap());
        let coeff = coeff_dist[random_index];
        coeff_dist.remove(random_index);
        let signed : bool = rng.gen();
        if (coeff > 0) && signed {
            coefficient_vec.push(-coeff);
        }
        else {
            coefficient_vec.push(coeff);
        }
    }
    let mut poly = R_q::new(coefficient_vec);
    poly
}


// Conjugation automorphism applied to a vector of polynomials
pub fn sigma_inv_vec(vec: &Vec<R_q>) -> Vec<R_q> {
    let mut res = vec![];
    for i in 0..vec.len() {
        res.push(sigma_inv(&vec[i]));
    }
    res
}



// The "Conjugation automorphism"
// This basically replaces all powers of X^n with X^{-n} and then reduces modulo X^d+1 (R's
// modulus) 
pub fn sigma_inv(a : &R_q) -> R_q {

    let poly_data_vec : Vec<Z_q> = a.data_vec();

    let mut new_coeff_vec : Vec<Z_q> = vec![Z_q::zero(); 64]; 

    // terms of degree 0 unchanged
    // degree 1-64, we simply transform X^n to -X^{64-n}.

    for deg in 0..poly_data_vec.len() {
        let coeff = poly_data_vec[deg];
        if deg == 0 {
            new_coeff_vec[0] = coeff;
        }
        else {
            let new_deg = (D as usize) - deg;
            new_coeff_vec[new_deg] = -coeff;
        }
    }
    let res_polynomial = R_q::new(new_coeff_vec);
    res_polynomial
}



pub fn multiply_poly_vec_ints(p_vec : &Vec<R_q>, ints: &Vec<Z_q>) -> Vec<R_q> {
    let mut p_vec_res : Vec<R_q> = vec![];

    for p in p_vec {
        p_vec_res.push(multiply_poly_ints(p, ints));
    }
    p_vec_res
}

pub fn multiply_poly_ints(p : &R_q, ints: &Vec<Z_q>) -> R_q {
    let mut p_res : R_q = R_q::new(vec![]);

    for coeff in ints {
        p_res = p_res + scale_polynomial(p, f32::from(*coeff));
    }
    p_res
}


// randomly samples n integers mod q, returns them as a vec
pub fn random_sample_Z_q(n : i128, q: i128) -> Vec<Z_q> {

    let mut rng = thread_rng();
    let dist = Uniform::from(0..q);

    let sample: Vec<Z_q> = Z_q::lift(&(0..n).map(|_| rng.sample(dist)).collect());
    sample
}

// scales a given polynomial by a scale factor (usually < 1)
pub fn scale_poly_vec(vec : &Vec<R_q>, s : f32) -> Vec<R_q> {
    let mut new_vec = vec![];
    for i in 0..vec.len() {
        new_vec.push(scale_polynomial(&vec[i], s));
    }
    new_vec
}


// scales a given polynomial by a scale factor (usually < 1)
pub fn scale_polynomial(p : &R_q, s : f32) -> R_q {
    //let constant_polynomial : R_q<f32> = R_q::new(vec![s]);
    let poly_vec : Vec<Z_q> = p.data_vec(); 
    let scaled_poly_vec : Vec<Z_q> = poly_vec.iter().map(|&x| (Z_q::from(((f32::from(x)) * s).floor()))).collect();
    return R_q::new(scaled_poly_vec);
}

pub fn scale_poly_vec_int(vec : &Vec<R_q>, s : &Z_q) -> Vec<R_q> {
    let mut new_vec = vec![];
    for i in 0..vec.len() {
        new_vec.push(scale_polynomial_int(&vec[i], s));
    }
    new_vec
}

pub fn scale_polynomial_rational(p : &R_q, a : &Z_q, b : &Z_q) -> R_q {
    let poly_vec : Vec<Z_q> = p.data_vec(); 
    let scaled_poly_vec : Vec<Z_q> = poly_vec.iter().map(|&x| ((x * *a) / *b)).collect();
    return R_q::new(scaled_poly_vec);
}


pub fn scale_polynomial_int(p : &R_q, s : &Z_q) -> R_q {
    let poly_vec : Vec<Z_q> = p.data_vec(); 
    let scaled_poly_vec : Vec<Z_q> = poly_vec.iter().map(|&x| x * *s).collect();
    return R_q::new(scaled_poly_vec);
}

pub fn vec_poly_norm_squared(vec: &Vec<R_q>) -> f64 {
    vec.iter()
        .map(|poly| poly_norm(poly))  // Compute the squared norm of each polynomial
        .sum()  // Sum the squared norms
}


// takes the 2-norm of a given polynomial (squared)
pub fn poly_norm(p : &R_q) -> f64 {
    let norm : Z_q = p.data_vec().iter().map(|&x| x*x).sum::<Z_q>();
    f64::from(norm)
}

// Does NOT square the norm
pub fn poly_norm_strict(p : &R_q) -> f64 {
    let strict_norm : f64 = f64::from(p.data_vec().iter().map(|&x| x*x).sum::<Z_q>()).sqrt();
    strict_norm
}


// Does NOT square the norm. Integer vec version
pub fn l2_norm(vec : &Vec<i128>) -> f64 {
    let sum_of_squares: i128 = vec.iter().map(|&x| x * x).sum();
    (sum_of_squares as f64).sqrt()
}

// TODO we're using a statistical estimate of the sup to accomplish this, not sure if that's valid.
// should be since there's basically no other way to feasibly do this (I think). Check this.
pub fn operator_norm(c : &R_q) -> f64 {
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

pub fn compute_total_norm(projection: &Array2<R_q>) -> f64 {
    let mut total_norm_squared: f64 = 0.0;
    for row in projection.outer_iter() {
        for poly in row.iter() {
            let norm = poly_norm(poly);
            total_norm_squared += norm;
        }
    }
    return f64::sqrt(total_norm_squared);
}

pub fn vec_inner_product(v1: &Vec<i128>, v2: &Vec<i128>) -> i128 {
    //println!("v1 length: {}", v1.len());
    //println!("v2 length: {}", v2.len());
    assert!(v1.len() == v2.len(), "inner product not defined on vectors of unequal length");
    let mut result : i128 = 0;
    for i in 0..v1.len() {
        let product = &v1[i] * &v2[i];
        result += product;
    }
    result
}
// TODO make these generic
pub fn vec_inner_product_Z_q(v1: &Vec<Z_q>, v2: &Vec<Z_q>) -> i128 {
    assert!(v1.len() == v2.len(), "inner product not defined on vectors of unequal length");
    let mut result : Z_q = Z_q::zero();
    for i in 0..v1.len() {
        let product = &v1[i] * &v2[i];
        result += product;
    }
    i128::from(result)
}


// elementwise multiplication of a vector of polynomials by a single polynomial. NOT an inner
// product.
pub fn poly_by_poly_vec(poly: &R_q, v1: &Vec<R_q>) -> Vec<R_q> {
    let mut v2 : Vec<R_q> = vec![];

    for i in 0..v1.len() {
        v2.push(poly * &v1[i]);
    }
    v2
}


pub fn add_poly_vec(v1: &Vec<R_q>, v2: &Vec<R_q>) -> Vec<R_q> {
    assert!(v1.len() == v2.len(), "summation not defined on vectors of unequal length");
    let mut v_res : Vec<R_q> = vec![];
    for i in 0..v1.len() {
        let new_p = &v1[i] + &v2[i];
        v_res.push(new_p);
    }
    v_res
}

// TODO genericize later
pub fn add_vecs(v1: &Vec<i128>, v2: &Vec<i128>) -> Vec<i128> {
    assert!(v1.len() == v2.len(), "summation not defined on vectors of unequal length");
    let mut v_res : Vec<i128> = vec![];
    for i in 0..v1.len() {
        let new_p = &v1[i] + &v2[i];
        v_res.push(new_p);
    }
    v_res
}



/* Rarely, we need to do elementwise addition of a single polynomial to an entire vector of
 * polynomials.. kind of strange, but how this is described in the protocol */
pub fn add_poly_vec_by_poly(v1: &Vec<R_q>, p: &R_q) -> Vec<R_q> {
    let mut v_res : Vec<R_q> = vec![];
    for i in 0..v1.len() {
        let new_p = &v1[i] + p;
        v_res.push(new_p);
    }
    v_res
}

pub fn gen_empty_poly_vec(n : usize) -> Vec<R_q> {
    let mut v_res : Vec<R_q> = vec![];
    for i in 0..n {
        v_res.push(R_q::new(vec![]));
    }

    v_res
}

// used specifically do decompose e.g., the entirety of t_i = t_i^(0) + ... + t_i^(t_1-1)b_1^(t_1-1)
pub fn decompose_polynomial_vec(vec : &Vec<R_q>, base : i128, exp: i128) -> Vec<Vec<R_q>> {

    // TODO KAPPA is hardcoded here for now. Fix later.
    let mut res = vec![vec![R_q::new(vec![]); KAPPA as usize]; exp as usize];

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
pub fn centered_rep(mut val : Z_q, b: i128) -> Z_q {
    //println!("centered representative computation! val: {}, b: {}", val, b);
    if val > b / 2 {
        //val -= Z_q::from(b);
        //TODO mod here is probably redundant
        return Z_q::new((b - i128::from(val)) % b);
    }
    //println!("val transformed to: {}", val);
    // no negative values are possible since this is Z_q, so we're done..
    val
}


pub fn decompose_polynomial(p : &R_q, base : i128, exp: i128) -> Vec<R_q> {

    let poly_data_vec : Vec<Z_q> = p.data_vec();

    // this takes the form of, e.g., g_{ij} = g_{ij}^(0) + ... + g_{ij}^{t_2-1}b_2^{t_2-1}
    // i.e., let's denote the whole polynomial a_j.
    // then each element of decomp_vec is the k-th component of each decomposed coefficient, i.e.,
    // a_j_k, found in the term a_j^(k)b^(k)
    let mut decomp_vec : Vec<R_q> = vec![R_q::new(vec![]); exp as usize];


    // has length deg(p)+1. Each element is an exp length vector..
    let mut decomposed_coeffs : Vec<R_q> = vec![R_q::new(vec![]); poly_data_vec.len()];

    // we decompose each coefficient a_j of the polynomial (which we can express as a base K 

    for deg in 0..poly_data_vec.len() {
        let mut poly_coeff = poly_data_vec[deg];
        //println!("deg: {}, coeff {}" , deg, poly_coeff);
        let mut decomposed_coefficient : Vec<Z_q> = vec![]; // coefficients of decomposed
                                                            // coefficient polynomial (soon to be)
        while poly_coeff != 0 {
            //println!("While loop!");
            let remainder = centered_rep(poly_coeff % base, base);
            decomposed_coefficient.push(remainder);
            poly_coeff = (poly_coeff - remainder) / base;
            //println!("remainder: {}. Poly coeff is now: {}", remainder, poly_coeff);
        }
        //println!("PHEW! done with that loop. Now decomposed_coefficient: {:?}", decomposed_coefficient);
        decomposed_coeffs[deg] = (R_q::new(decomposed_coefficient));
    }
    //println!("DONE!");

    // Next, we transform the "decomposed_coeffs" into the actual decomp_vec by taking the k-th
    // entries from each coefficient and putting them together..
    for k in 0..exp {
        let mut a_j_k : R_q = R_q::zero();
        //println!("Decomposed coeffs len: {}", decomposed_coeffs.len());
        for i in 0..decomposed_coeffs.len() {
            //println!("decomposed coeff i: {}, k={}, polynomial: {}", i,k, &decomposed_coeffs[i]);
            a_j_k = a_j_k + (&decomposed_coeffs[i]).get_term_of_deg(k as usize);
        }
        decomp_vec[k as usize] = a_j_k;
    }
    decomp_vec
}

// Used to convert a vector of polynomials \vec{\b{s}}_i to another vector
// which is a concatenation of all the coefficients 
pub fn witness_coeff_concat(vec: &Vec<R_q>) -> Vec<Z_q> {
    // NOTE: we assume that the polynomials in vec will be of degree D... 
    // But we want to return a vec of degree N*D (where N is the assumed length of vec)..
    // So we recognize that for a polynomial ring mod X^D+1, X^d congruent to -1... so we can
    // rewrite as such.
    let mut coeffs : Vec<Z_q> = vec![];
    for j in 0..vec.len() {
        let poly : &R_q = &vec[j];
        let poly_vec_data : Vec<Z_q> = poly.data_vec();
        
        // TODO can we simplify this to poly_vec_data.len()? I don't think so, since data might
        // store fewer coeffs
        for deg in 0..(D as usize) {
            if deg < poly_vec_data.len() {
                coeffs.push(poly_vec_data[deg]);
            }
            else {
                coeffs.push(Z_q::zero());
            }
        }
    }
    coeffs
}

// Basically does the opposite of the above function
// Takes a vec of straight ints of dimension N*D for some int N, and returns
// N length vec of reconstructed R_qs using the coefficients
pub fn concat_coeff_reduction(vec: &Vec<Z_q>) -> Vec<R_q> {
    assert!((vec.len() % (D as usize)) == 0, "dimension of vec is wrong");
    let n = vec.len() / (D as usize);
    let mut res = vec![];
    for i in (0..(n*(D as usize))).step_by(D as usize) {
        let new_vec = vec[i..(i+D as usize)].to_vec();         
        let poly = R_q::new(new_vec);
        res.push(poly);
    }
    res
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
pub fn polynomial_vec_inner_product(v1: &[R_q], v2: &[R_q]) -> R_q {
    assert!(v1.len() == v2.len(), "inner product not defined on vectors of unequal length. v1 length: {}, v2 length: {}", v1.len(), v2.len());
    //println!("\n\nNew Inner product!");
    let mut result = R_q::new(vec![]);
    for i in 0..v1.len() {
        let product = &v1[i] * &v2[i];
        /*
        if (product == R_q::new(vec![])) { 
            println!("PRODUCT ZERO!");
            println!("v1: {}", &v1[i]);
            println!("v2: {}", &v2[i]);
        }
        else {
            println!("NOT ZERO!!!");
        }
        */
        result = result + product;
    }
    result
}

pub fn matmul(m1: &Array2<i128>, m2: &Array2<i128>) -> Array2<i128> {
    assert!(m1.shape()[1] == m2.shape()[0], "matrix product not defined on matrices which do not have dimension (m,n) x (n,k), for some m,n,k");
    //println!("m1 shape: {:?} \n m2 shape: {:?}", m1.shape(), m2.shape());
    let m = m1.shape()[0];
    let k = m2.shape()[1];

    let mut result = Array2::from_elem(Ix2(m, k), 0); 
    for i in 0..m {
        for j in 0..k {
            let v1 = m1.row(i).to_vec();
            let v2 = m2.column(j).to_vec();
            result[[i,j]] = vec_inner_product(&v1, &v2);
        }
    }
    result
}

// compute the matrix product between two polynomial MATRICES (Array2, slightly different)
pub fn polynomial_matrix_product(m1: &Array2<R_q>, m2: &Array2<R_q>) -> Array2<R_q> {
    //println!("m1 shape: {:?} \n m2 shape: {:?}", m1.shape(), m2.shape());
    assert!(m1.shape()[1] == m2.shape()[0], "matrix product not defined on matrices which do not have dimension (m,n) x (n,k), for some m,n,k");
    let m = m1.shape()[0];
    let k = m2.shape()[1];

    let mut result = Array2::from_elem(Ix2(m, k), R_q::new(vec![])); 
    for i in 0..m {
        for j in 0..k {
            let v1 = m1.row(i).to_vec();
            let v2 = m2.column(j).to_vec();
            result[[i,j]] = polynomial_vec_inner_product(&v1, &v2);
        }
    }
    result
}

pub fn sample_jl_projection_gen() -> Array2<i128> {

    let between = Uniform::from(-1..=1);

    let mut rng = rand::thread_rng();

    let mut Pi_i : Array2<i128> = Array2::zeros((256, N*(D as usize)));

    for ((i, j), value) in Pi_i.indexed_iter_mut() {
        *value = between.sample(&mut rng);
    }


    // TODO don't entirely understand the functionality of this line.. but seems to work.
    Pi_i 
}

pub fn jl_project_gen(vec: &Vec<i128>) -> Vec<i128> {
    //let between = Uniform::from(-1..=1);
    let mut rng = rand::thread_rng();
    let choices = [-1,0,1];
    let weights = [0.25, 0.5, 0.25];
    let dist = WeightedIndex::new(&weights).unwrap();



    // we get random matrices in {-1,0,1}
    let mut projection : Vec<i128> = vec![0 ; 256];
    //for i in 0..R {
        //let Pi_i = sample_jl_projection_gen();
        let mut Pi_i : Array2<i128> = Array2::zeros((256, N*(D as usize)));
        for ((i, j), value) in Pi_i.indexed_iter_mut() {
            *value = choices[dist.sample(&mut rng)] as i128;
        }
        //println!("Got Pi_i for i={}",i);
        let s_i_coeffs : Array2<i128> = vec_to_column_array(vec);
        println!("coeffs: {:?}", s_i_coeffs.column(0).to_vec());
        println!("Pi_i column 0: {:?}", Pi_i.column(0).to_vec());
        //println!("Got s_i coeffs");
        // NOTE: for reference, this is a 256x(ND) multiplied by an (ND)x1, giving a 256x1
        // which we turn into a vec
        let product = matmul(&Pi_i, &s_i_coeffs).column(0).to_vec();
        println!("JL PRODUCT! : {:?}", product);
        //println!("computed product");
        projection = add_vecs(&projection, &product);
    //}
    projection
}


