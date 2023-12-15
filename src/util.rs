use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

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


pub fn generate_polynomial_picky(q : i64, d : i64, coeff_dist : vec<i64>) -> Polynomial<i64> {
    assert!(coeff_dist.len() == d, "Must have one coefficient for each degree of polynomial");
    let mut rng = rand::thread_rng();
    let mut coefficient_vec: Vec<i64> = Vec::new();
    for degree in 0..d {
        let random_index = coeff_dist.as_slice().choose(&mut rng).unwrap();
        let coeff = coeff_dist[*random_index];
        coeff_dist.remove(*random_index);
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

// The Legendary "Conjugation automorphism"
// Computes the ring inverse of a polynomial element of R_q
pub fn sigma_inv(a : Polynomial) -> Polynomial {
    // TODO actually implement
    a
}



pub fn multiply_poly_vec_ints(p_vec : Vec<Polynomial<i64>>, ints: Vec<i64>) -> Vec<Polynomial<i64>> {
    let p_vec_res : Vec<Polynomial<i64>> = vec![];

    for p in &p_vec {
        p_vec_res.push(multiply_poly_ints(p, ints));
    }
    p_vec_res
}

pub fn multiply_poly_ints(p : Polynomial<i64>, ints: Vec<i64>) -> Polynomial<i64> {
    let p_res : Polynomial<i64> = Polynomial::new(vec![]);

    for coeff in &ints {
        p_res += scale_polynomial(p, coeff as f32);
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
pub fn scale_polynomial(p : Polynomial<i64>, s : f32) -> Polynomial<i64> {
    //let constant_polynomial : Polynomial<f32> = Polynomial::new(vec![s]);
    let poly_vec : Vec<i64> = p.data().to_vec(); 
    let scaled_poly_vec : Vec<i64> = poly_vec.iter().map(|&x| (((x as f32) * s).floor() as i64)).collect();
    return Polynomial::new(scaled_poly_vec);
}

// takes the 2-norm of a given polynomial (squared)
pub fn poly_norm(p : Polynomial<i64>) -> f64 {
    let norm : i64 = p.data().to_vec().iter().map(|&x| x*x).sum::<i64>();
    norm as f64
}


// TODO this is terrible. Zero idea how to do this one at the moment.
// But I'm certain it can be done. Ask around, think about it. Do it later.
pub fn operator_norm(p : Polynomial<i64>) -> f64 {
    10.0
}





pub fn compute_total_norm(projection: Array2<Polynomial<i64>>) -> f64 {

    let mut total_norm_squared: f64 = 0.0;


    for row in projection.outer_iter() {
        for poly in row.iter() {
            let norm = poly_norm(poly.clone());
            total_norm_squared += norm;
        }
    }
    return f64::sqrt(total_norm_squared);
}

pub fn vec_inner_product(v1: Vec<i64>, v2: Vec<i64>) -> Vec<i64> {
    assert!(v1.len() == v2.len(), "inner product not defined on vectors of unequal length");
    let mut result = vec![];
    for i in 0..v1.len() {
        let product = &v1[i] * &v2[i];
        result = result + product;
    }
    result
}


// elementwise multiplication of a vector of polynomials by a single polynomial. NOT an inner
// product.
pub fn poly_by_poly_vec(poly: Polynomial<i64>, v1: Vec<Polynomial<i64>>) -> Vec<Polynomial<i64>> {
    let mut v2 : Vec<Polynomial<i64>> = vec![];

    for i in 0..v1.len() {
        v2.push(poly * v1[i]);
    }
    v2
}






// compute the dot product between two lists of polynomials
pub fn polynomial_vec_inner_product(v1: Vec<Polynomial<i64>>, v2: Vec<Polynomial<i64>>) -> Polynomial<i64> {
    assert!(v1.len() == v2.len(), "inner product not defined on vectors of unequal length");
    let mut result = Polynomial::new(vec![]);
    for i in 0..v1.len() {
        let product = &v1[i] * &v2[i];
        result = result + product;
    }
    result
}

// compute the matrix product between two polynomial MATRICES (Array2, slightly different)
pub fn polynomial_matrix_product(m1: Array2<Polynomial<i64>>, m2: Array2<Polynomial<i64>>) -> Array2<Polynomial<i64>> {
    assert!(m1.shape()[1] == m2.shape()[0], "matrix product not defined on matrices which do not have dimension (m,n) x (n,k), for some m,n,k");
    let m = m1.shape()[0];
    let k = m2.shape()[1];

    let mut result = Array2::from_elem(Ix2(m, k), Polynomial::new(vec![])); 
    for i in 0..m {
        for j in 0..k {
            let v1 = m1.row(i).to_vec();
            let v2 = m1.column(j).to_vec();
            result[[i,j]] = polynomial_vec_inner_product(v1, v2);
        }
    }
    result
}
