use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

// generates an element of Z_q[X]/(X^d+1)
// TODO: more sophisticated / cleaner random polynomial sampling
pub fn generate_polynomial(q : i64, d : i64) -> Polynomial<i64> {
    let max_degree = d-1;
    let mut rng = rand::thread_rng();
    //let num_terms = rng.gen_range(1, max_degree + 2);
    let mut coefficient_vec: Vec<i64> = Vec::new();
    for degree in 0..max_degree {
        coefficient_vec.push(rng.gen_range(0..q));
    }
    let mut poly = Polynomial::new(coefficient_vec);
    poly
}

/*
pub fn jl_projection(w: &[i64]) -> Vec<i64> {

    // proving knowledge of a long vector w in Z^d without revealing it

    // Let verifier sample a random linear map Pi: Z^d \to Z^256
    // entries of Pi are independent and equal to -1,0,1 with probabilities 1/4, 1/2, 1/4

    // sample random map Pi

    let d = w.len();
    let mut projection_vec = vec![0; 256];

    for i in 0..255 {
        let mut random_buffer = vec![0;d];
        for j in 0..(d-1) {
            if rand::random() {
                if rand::random() {
                    random_buffer[i] = 1;
                }
                else {
                    random_buffer[i] = -1;
                }
            }
            else {
                random_buffer[i] = 0;
            }
        }
        let v1 = DVector::from_vec(random_buffer);
        let v2 = DVector::from_row_slice(w);
        projection_vec[i] = v1.dot(&v2);
    }
    return projection_vec;
}
*/


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
pub fn poly_norm(p : Polynomial<i64>) -> i64 {
    let norm : i64 = p.data().to_vec().iter().map(|&x| x*x).sum::<i64>();
    norm
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
