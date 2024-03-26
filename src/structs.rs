use ndarray::Array2;
use std::collections::HashMap;
use std::mem;

use crate::algebraic::*;
use crate::constants::*;
use crate::util::*;

pub struct CRS {
    pub a_mat: Array2<Rq>,
    pub b_mat: HashMap<(usize, usize), Array2<Rq>>,
    pub c_mat: HashMap<(usize, usize, usize), Array2<Rq>>,
    pub d_mat: HashMap<(usize, usize, usize), Array2<Rq>>,
}

impl CRS {
    // TODO random uniform generation for now... perhaps this warrants something more sophisticated
    // later.
    pub fn new(constants: &RuntimeConstants) -> Self {
        if is_verbose() {
            println!("Generating random matrix A");
        }
        let a_mat = generate_random_matrix(constants.KAPPA, constants.N, *Q, D);
        /*
        for i in 0..KAPPA {
            for j in 0..N {
                println!("A val at i={}, j={}: {}", i, j, &A_mat[[i,j]]);
            }
        }
        */
        let mut b_mat: HashMap<(usize, usize), Array2<Rq>> = HashMap::new();
        let mut c_mat: HashMap<(usize, usize, usize), Array2<Rq>> = HashMap::new();
        let mut d_mat: HashMap<(usize, usize, usize), Array2<Rq>> = HashMap::new();

        if is_verbose() {
            println!("Generating B matrices");
        }

        // TODO maybe put this all in one loop to clean it up, but then again maybe not.
        // TODO it's technically 1 \leq i \leq j \leq R, i.e., 1..(R+1), but... we're not doing
        // that, at least for now.
        for i in 0..constants.R {
            for k in 0..(constants.T_1 as usize) {
                let b_ik = generate_random_matrix(constants.KAPPA_1, constants.KAPPA, *Q, D);
                let index = (i, k);
                b_mat.insert(index, b_ik);
            }
        }

        if is_verbose() {
            println!("Generating C matrices");
        }

        for i in 0..constants.R {
            for j in i..constants.R {
                for k in 0..(constants.T_2 as usize) {
                    let c_ijk = generate_random_matrix(constants.KAPPA_2, 1, *Q, D);
                    let index = (i, j, k);
                    c_mat.insert(index, c_ijk);
                }
            }
        }

        if is_verbose() {
            println!("Generating D matrices");
        }

        for i in 0..constants.R {
            for j in i..constants.R {
                for k in 0..(constants.T_1 as usize) {
                    let d_ijk = generate_random_matrix(constants.KAPPA_2, 1, *Q, D);
                    let index = (i, j, k);
                    d_mat.insert(index, d_ijk);
                }
            }
        }
        CRS {
            a_mat,
            b_mat,
            c_mat,
            d_mat,
        }
    }
}

pub struct Transcript {
    // fields (see protocol)
    pub u_1: Vec<Rq>,
    pub pi_i_all: Vec<Array2<Zq>>,
    pub projection: Vec<Zq>,
    pub psi: Vec<Vec<Zq>>,   // note: This contains all ceil(128/log(q)) psi_k
    pub omega: Vec<Vec<Zq>>, // note: This contains all ceil(128/log(q)) omega_k
    pub b_prime_prime: Vec<Rq>,
    pub alpha: Vec<Rq>,
    pub beta: Vec<Rq>,
    pub u_2: Vec<Rq>,
    pub c: Vec<Rq>,
    pub z: Vec<Rq>,
    pub t_i_all: Vec<Vec<Rq>>,
    pub g_mat: Array2<Rq>,
    pub h_mat: Array2<Rq>,
}

impl Transcript {
    pub fn size_in_bytes(&self) -> usize {
        let mut size = 0;
        size += mem::size_of_val(&self.u_1);
        size += mem::size_of_val(&self.pi_i_all);
        size += mem::size_of_val(&self.projection);
        size += mem::size_of_val(&self.psi);
        size += mem::size_of_val(&self.omega);
        size += mem::size_of_val(&self.b_prime_prime);
        size += mem::size_of_val(&self.alpha);
        size += mem::size_of_val(&self.beta);
        size += mem::size_of_val(&self.u_2);
        size += mem::size_of_val(&self.c);
        size += mem::size_of_val(&self.z);
        size += mem::size_of_val(&self.t_i_all);
        size += mem::size_of_val(&self.g_mat);
        size += mem::size_of_val(&self.h_mat);
        size
    }
}

pub struct State {
    // This is how we store the families of functions F and F'...
    // For now, we just put them in a vec in the State.. the ith index
    // stores the corresponding matrix for the ith function in the principal relation.

    // TODO: for now, we'll also just keep K=1, L=1, and we'll have F' consist of the one function
    // in F.

    // Contents of F:
    pub phi_k: Vec<Array2<Rq>>,
    pub a_k: Vec<Array2<Rq>>,
    pub b_k: Vec<Rq>,

    // Contents of F' (constant term of the solution must be zero):
    pub phi_prime_k: Vec<Array2<Rq>>,
    pub a_prime_k: Vec<Array2<Rq>>,
    pub b_prime_k: Vec<Zq>,
}

impl State {
    fn gen_f(witness: &Array2<Rq>, constants: &RuntimeConstants) -> (Array2<Rq>, Array2<Rq>, Rq) {
        // random generation of the polynomial matrix a_constraints (NOT TO BE CONFUSED WITH CRS matrix A!)
        // TODO maybe we don't need "constants" here, we can just fetch R and N from the witness
        // since that's really all we need?
        let mut a_constraints = Array2::from_elem((constants.R, constants.R), Rq::new(vec![]));
        for i in 0..constants.R {
            for j in 0..constants.R {
                let a_ij = generate_polynomial(*Q, D);
                if a_constraints[[i, j]] == Rq::new(vec![]) {
                    a_constraints[[i, j]] = a_ij.clone();
                    a_constraints[[j, i]] = a_ij;
                }
            }
        }
        if is_verbose() {
            println!("Generated a_constraints!");
        }

        // random generation of random polynomial matrix phi
        let mut phi = Array2::from_elem((constants.N, constants.R), Rq::new(vec![]));
        for i in 0..constants.R {
            for j in 0..constants.N {
                let phi_ji = generate_polynomial(*Q, D);
                phi[[j, i]] = phi_ji;
            }
        }

        // Next, we need to compute 'b' such that this relation is equal to zero
        // So we compute the values of the relation and work backwards
        let mut a_product: Rq = Rq::new(vec![]);
        for i in 0..constants.R {
            for j in 0..constants.R {
                let inner_prod = polynomial_vec_inner_product(
                    &witness.column(i).to_vec(),
                    &witness.column(j).to_vec(),
                );
                let prod = a_constraints[[i, j]].clone() * inner_prod;
                a_product = a_product.clone() + prod;
            }
        }

        let mut phi_product: Rq = Rq::new(vec![]);
        for i in 0..constants.R {
            let inner_prod =
                polynomial_vec_inner_product(&phi.column(i).to_vec(), &witness.column(i).to_vec());
            phi_product = phi_product.clone() + inner_prod;
        }

        let b: Rq = &a_product + &phi_product;
        //println!("Generated b!");
        //println!("{}\n\n", b);

        //println!("A product: {}\n", a_product);
        //println!("phi product: {}", phi_product);

        (phi, a_constraints, b)
    }

    pub fn new(witness: &Array2<Rq>, constants: &RuntimeConstants) -> Self {
        let mut phi_k: Vec<Array2<Rq>> = vec![];
        let mut a_k: Vec<Array2<Rq>> = vec![];
        let mut b_k: Vec<Rq> = vec![];

        let mut phi_prime_k: Vec<Array2<Rq>> = vec![];
        let mut a_prime_k: Vec<Array2<Rq>> = vec![];
        let mut b_prime_k: Vec<Zq> = vec![];

        for _k in 0..K {
            let (phi, a, b) = Self::gen_f(witness, constants);

            // TODO obviously cloning is not great here, but references are an even worse option
            // because no immediate way to do this without creating dangling pointer / not
            // compiling.. so we will potentially refactor later once minimally working.
            phi_k.push(phi.clone());
            a_k.push(a.clone());
            b_k.push(b.clone());

            phi_prime_k.push(phi);
            a_prime_k.push(a);
            b_prime_k.push(b.eval(Zq::from(0)));
        }
        for _l in 0..L {
            // generate corresponding f' functions.
            // TODO actually do something here, maybe
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
