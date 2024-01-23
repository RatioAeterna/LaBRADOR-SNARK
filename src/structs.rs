use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;
use std::collections::HashMap;

use crate::util::*;
use crate::constants::*;


pub struct CRS {
    // TODO should these fields be public? Perhaps not. Will fix later. 
    pub A_mat : Array2<Polynomial<i64>>,
    pub B_mat : HashMap<(usize, usize), Array2<Polynomial<i64>>>,
    pub C_mat : HashMap<(usize, usize, usize), Array2<Polynomial<i64>>>,
    pub D_mat : HashMap<(usize, usize, usize), Array2<Polynomial<i64>>>,
}


impl CRS {
    // TODO random uniform generation for now... perhaps this warrants something more sophisticated
    // later.
    pub fn new() -> Self {
        let A_mat = generate_random_matrix(KAPPA, N as i64, Q, D);
        let mut B_mat : HashMap<(usize, usize), Array2<Polynomial<i64>>> = HashMap::new();
        let mut C_mat : HashMap<(usize, usize, usize), Array2<Polynomial<i64>>> = HashMap::new();
        let mut D_mat : HashMap<(usize, usize, usize), Array2<Polynomial<i64>>> = HashMap::new();

        // TODO maybe put this all in one loop to clean it up, but then again maybe not.
        for i in 1..(R+1) {
            for k in 0..(*T_1 as usize) {
                let B_ik = generate_random_matrix(KAPPA_1, N as i64, Q, D);
                let index = (i,k);
                B_mat.insert(index, B_ik);
            }
        }

        for i in 1..(R+1) {
            for j in i..(R+1) {
                for k in 0..(*T_2 as usize) {
                    let C_ijk = generate_random_matrix(KAPPA_2, 1, Q, D);
                    let index = (i,j,k);
                    C_mat.insert(index, C_ijk);
                }
            }
        }

        for i in 1..(R+1) {
            for j in i..(R+1) {
                for k in 0..(*T_1 as usize) {
                    let D_ijk = generate_random_matrix(KAPPA_2, 1, Q, D);
                    let index = (i,j,k);
                    D_mat.insert(index, D_ijk);
                }
            }
        }
        CRS {A_mat, B_mat, C_mat, D_mat, }
    }
}

pub struct Transcript {
    // fields (see protocol)
    pub u_1 : Vec<Polynomial<i64>>,
    pub projection : Vec<i64>,
    pub psi : Vec<Vec<i64>>, // note: This contains all ceil(128/log(q)) psi_k
    pub omega : Vec<Vec<i64>>, // note: This contains all ceil(128/log(q)) omega_k
    pub alpha : Vec<Polynomial<i64>>,
    pub beta : Vec<Polynomial<i64>>,
    pub u_2 : Vec<Polynomial<i64>>,
    pub c : Vec<Polynomial<i64>>,
    pub z : Vec<Polynomial<i64>>,
    pub Gij : Array2<Polynomial<i64>>,
    pub Hij : Array2<Polynomial<i64>>,
}


pub struct State {
    // This is how we store the families of functions F and F'... 
    // For now, we just put them in a vec in the State.. the ith index
    // stores the corresponding matrix for the ith function in the principal relation.

    // TODO: for now, we'll also just keep K=1, L=1, and we'll have F' consist of the one function
    // in F.

    // Contents of F:
    pub phi_k : Vec<Array2<Polynomial<i64>>>,
    pub a_k : Vec<Array2<Polynomial<i64>>>,
    pub b_k : Vec<Polynomial<i64>>,

    // Contents of F' (constant term of the solution must be zero):
    pub phi_prime_k : Vec<Array2<Polynomial<i64>>>,
    pub a_prime_k : Vec<Array2<Polynomial<i64>>>,
    // TODO these should be evaluations at ZERO, i.e., i64 constant coefficients, not polynomials
    pub b_prime_k : Vec<Polynomial<i64>>,
}

impl State {

    fn gen_f(S: &Array2<Polynomial<i64>>) -> (Array2<Polynomial<i64>>, Array2<Polynomial<i64>>, Polynomial<i64>) {
        // random generation of the polynomial matrix Aij (NOT TO BE CONFUSED WITH CRS matrix A!)
        let mut Aij = Array2::from_elem((R,R), Polynomial::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let a_ij = generate_polynomial(Q,D);
                if Aij[[i,j]] == Polynomial::new(vec![]) {
                    Aij[[i,j]] = a_ij.clone();
                    Aij[[j,i]] = a_ij;
                }
            }
        }
        println!("Generated Aij!");
        for row in Aij.rows() {
            for poly in row {
                println!("{}", poly.pretty("x"));
            }
        }

        // random generation of random polynomial matrix Phi
        let mut Phi = Array2::from_elem((N,R), Polynomial::new(vec![])); 
        for i in 0..N {
            for j in 0..R {
                let phi_ij = generate_polynomial(Q,D);
                Phi[[i,j]] = phi_ij;
            }
        }

        // Next, we need to compute 'b' such that this relation is equal to zero
        // So we compute the values of the relation and work backwards
        let mut a_product : Polynomial<i64> = Polynomial::new(vec![]);
        for i in 0..R {
            for j in 0..R {
                let inner_prod = polynomial_vec_inner_product(&S.column(i).to_vec(), &S.column(j).to_vec());
                let prod = Aij[[i,j]].clone() * inner_prod;
                a_product = a_product.clone() + prod;
            }
        }

        let mut phi_product : Polynomial<i64> = Polynomial::new(vec![]);
        for i in 0..R {
            let inner_prod = polynomial_vec_inner_product(&Phi.column(i).to_vec(), &S.column(i).to_vec());
            phi_product = phi_product.clone() + inner_prod;
        }

        let b : Polynomial<i64> = a_product + phi_product;
        println!("Generated b!");
        println!("{}", b.pretty("x"));

        (Phi, Aij, b)
    }




    pub fn new(S: &Array2<Polynomial<i64>>) -> Self {
        let mut phi_k : Vec<Array2<Polynomial<i64>>> = vec![];
        let mut a_k : Vec<Array2<Polynomial<i64>>> = vec![];
        let mut b_k : Vec<Polynomial<i64>> = vec![];

        let mut phi_prime_k : Vec<Array2<Polynomial<i64>>> = vec![];
        let mut a_prime_k : Vec<Array2<Polynomial<i64>>> = vec![];
        let mut b_prime_k : Vec<Polynomial<i64>> = vec![];

        for k in 0..K {
            let phi : Array2<Polynomial<i64>>;
            let a : Array2<Polynomial<i64>>;
            let b : Polynomial<i64>;

            let (phi, a, b) = Self::gen_f(S);

            // TODO obviously cloning is not great here, but references are an even worse option
            // because no immediate way to do this without creating dangling pointer / not
            // compiling.. so we will potentially refactor later once minimally working.
            phi_k.push(phi.clone());
            a_k.push(a.clone());
            b_k.push(b.clone());

            phi_prime_k.push(phi);
            a_prime_k.push(a);
            b_prime_k.push(b);
        }
        for l in 0..L {
            // generate corresponding f' functions.
        }
        State {
            phi_k,
            a_k,
            b_k,
            phi_prime_k,
            a_prime_k,
            b_prime_k,
        }
    }
}



