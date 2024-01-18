use ndarray::{Array2, Ix2, concatenate};
use polynomial::Polynomial;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use rand::distributions::Uniform;

use crate::util::*;
use crate::constants::*;


pub struct CRS {
    // TODO should these fields be public? Perhaps not. Will fix later. 
    pub A : Array2<Polynomial<i64>>,
    pub B : HashMap<Array2<Polynomial<i64>>,
    pub C : HashMap<Array2<Polynomial<i64>>,
    pub D : HashMap<Array2<Polynomial<i64>>,
}


impl CRS {
    // TODO random uniform generation for now... perhaps this warrants something more sophisticated
    // later.
    pub fn new() -> Self {
        let A = generate_random_matrix(KAPPA, N as i64, Q, D);
        let mut B : HashMap<Array2<Polynomial<i64>> = HashMap::new();
        let mut C : HashMap<Array2<Polynomial<i64>> = HashMap::new();
        let mut D : HashMap<Array2<Polynomial<i64>> = HashMap::new();

        // TODO maybe put this all in one loop to clean it up, but then again maybe not.
        for i in 1..(R+1) {
            for k in 0..T_1 {
                B_ik = generate_random_matrix(KAPPA_1, N as i64, Q, D);
                let index = (i,k);
                B.insert(index, B_ik);
            }
        }

        for i in 1..(R+1) {
            for j in i..(R+1) {
                for k in 0..T_2 {
                    C_ijk = generate_random_matrix(KAPPA_2, 1, Q, D);
                    let index = (i,j,k);
                    C.insert(index, C_ijk);
                }
            }
        }

        for i in 1..(R+1) {
            for j in i..(R+1) {
                for k in 0..T_1 {
                    D_ijk = generate_random_matrix(KAPPA_2, 1, Q, D);
                    let index = (i,j,k);
                    D.insert(index, D_ijk);
                }
            }
        }
        CRS {A, B, C, D, }
    }
}

pub struct Transcript {
    // fields (see protocol)
    pub projection : Array2<Polynomial<i64>>,
    pub psi : Vec<Vec<i64>> // note: This contains all ceil(128/log(q)) psi_k
    pub omega : Vec<Vec<i64>> // note: This contains all ceil(128/log(q)) omega_k
    pub alpha : Vec<Polynomial<i64>>,
    pub beta : Vec<Polynomial<i64>>,
    pub u_2 : Vec<Polynomial<i64>>,
    pub c : Vec<Polynomial<i64>>,
    pub z : Vec<Polynomial<i64>>,
    pub gij : Array2<Polynomial<i64>>,
    pub hij : Array2<Polynomial<i64>>,
}


pub struct State {
    // This is how we store the families of functions F and F'... 
    // For now, we just put them in a vec in the State.. the ith index
    // stores the corresponding matrix for the ith function in the principal relation.

    // TODO: for now, we'll also just keep K=1, L=1, and we'll have F' consist of the one function
    // in F.

    // Contents of F:
    phi_k : Vec<Array2<Polynomial<i64>>,
    a_k : Vec<Array2<Polynomial<i64>>,
    b_k : Vec<Polynomial<i64>,

    // Contents of F' (constant term of the solution must be zero):
    phi_prime_k : Vec<Array2<Polynomial<i64>>,
    a_prime_k : Vec<Array2<Polynomial<i64>>,
    b_prime_k : Vec<Polynomial<i64>,
}

impl State {

    fn gen_f(S: Array2<Polynomial<i64>>) -> (Array2<Polynomial<i64>>, Array2<Polynomial<i64>>, Polynomial<i64>) {
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
                let inner_prod = polynomial_vec_inner_product(S.column(i).to_vec(), S.column(j).to_vec());
                let prod = Aij[[i,j]].clone() * inner_prod;
                a_product = a_product.clone() + prod;
            }
        }

        let mut phi_product : Polynomial<i64> = Polynomial::new(vec![]);
        for i in 0..R {
            let inner_prod = polynomial_vec_inner_product(Phi.column(i).to_vec(), S.column(i).to_vec());
            phi_product = phi_product.clone() + inner_prod;
        }

        let b : Polynomial<i64> = a_product + phi_product;
        println!("Generated b!");
        println!("{}", b.pretty("x"));

        (Phi, Aij, b)
    }




    pub fn new(S: Array2<Polynomial<i64>>) -> Self {
        let mut phi_k : Vec<Array2<Polynomial<i64>> = vec![];
        let mut a_k : Vec<Array2<Polynomial<i64>> = vec![];
        let mut b_k : Vec<Polynomial<i64> = vec![];

        for k in 0..K {
            let phi : Array2<Polynomial<i64>>;
            let a : Array2<Polynomial<i64>>;
            let b : Polynomial<i64>;

            let (phi, a, b) = gen_f(S)

            phi_k.push(phi);
            a_k.push(a);
            b_k.push(b);

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



