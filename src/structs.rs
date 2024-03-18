use ndarray::Array2;
use std::collections::HashMap;
use std::mem;

use crate::algebraic::*;
use crate::util::*;
use crate::constants::*;


pub struct CRS {
    // TODO should these fields be public? Perhaps not. Will fix later. 
    pub A_mat : Array2<Rq>,
    pub B_mat : HashMap<(usize, usize), Array2<Rq>>,
    pub C_mat : HashMap<(usize, usize, usize), Array2<Rq>>,
    pub D_mat : HashMap<(usize, usize, usize), Array2<Rq>>,
}


impl CRS {
    // TODO random uniform generation for now... perhaps this warrants something more sophisticated
    // later.
    pub fn new() -> Self {
        println!("Generating random matrix A");
        let A_mat = generate_random_matrix(KAPPA, N, *Q, D);
        /*
        for i in 0..KAPPA {
            for j in 0..N {
                println!("A val at i={}, j={}: {}", i, j, &A_mat[[i,j]]);
            }
        }
        */
        let mut B_mat : HashMap<(usize, usize), Array2<Rq>> = HashMap::new();
        let mut C_mat : HashMap<(usize, usize, usize), Array2<Rq>> = HashMap::new();
        let mut D_mat : HashMap<(usize, usize, usize), Array2<Rq>> = HashMap::new();

        println!("Generating B matrices");
        

        // TODO maybe put this all in one loop to clean it up, but then again maybe not.
        // TODO it's technically 1 \leq i \leq j \leq R, i.e., 1..(R+1), but... we're not doing
        // that, at least for now.
        for i in 0..R {
            for k in 0..(*T_1 as usize) {
                let B_ik = generate_random_matrix(KAPPA_1, KAPPA, *Q, D);
                let index = (i,k);
                B_mat.insert(index, B_ik);
            }
        }

        println!("Generating C matrices");

        for i in 0..R {
            for j in i..R {
                for k in 0..(*T_2 as usize) {
                    let C_ijk = generate_random_matrix(KAPPA_2, 1, *Q, D);
                    let index = (i,j,k);
                    C_mat.insert(index, C_ijk);
                }
            }
        }

        println!("Generating D matrices");

        for i in 0..R {
            for j in i..R {
                for k in 0..(*T_1 as usize) {
                    let D_ijk = generate_random_matrix(KAPPA_2, 1, *Q, D);
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
    pub u_1 : Vec<Rq>,
    pub Pi_i_all : Vec<Array2<Zq>>, 
    pub projection : Vec<Zq>,
    pub psi : Vec<Vec<Zq>>, // note: This contains all ceil(128/log(q)) psi_k
    pub omega : Vec<Vec<Zq>>, // note: This contains all ceil(128/log(q)) omega_k
    pub b_prime_prime: Vec<Rq>,
    pub alpha : Vec<Rq>,
    pub beta : Vec<Rq>,
    pub u_2 : Vec<Rq>,
    pub c : Vec<Rq>,
    pub z : Vec<Rq>,
    pub t_i_all : Vec<Vec<Rq>>, 
    pub Gij : Array2<Rq>,
    pub Hij : Array2<Rq>,
}


impl Transcript {

    pub fn size_in_bytes(&self) -> usize {
        let mut size = 0;
        size += mem::size_of_val(&self.u_1);
        size += mem::size_of_val(&self.Pi_i_all);
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
        size += mem::size_of_val(&self.Gij);
        size += mem::size_of_val(&self.Hij);
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
    pub phi_k : Vec<Array2<Rq>>,
    pub a_k : Vec<Array2<Rq>>,
    pub b_k : Vec<Rq>,

    // Contents of F' (constant term of the solution must be zero):
    pub phi_prime_k : Vec<Array2<Rq>>,
    pub a_prime_k : Vec<Array2<Rq>>,
    pub b_prime_k : Vec<Zq>,
}

impl State {

    fn gen_f(S: &Array2<Rq>) -> (Array2<Rq>, Array2<Rq>, Rq) {
        // random generation of the polynomial matrix Aij (NOT TO BE CONFUSED WITH CRS matrix A!)
        let mut Aij = Array2::from_elem((R,R), Rq::new(vec![])); 
        for i in 0..R {
            for j in 0..R {
                let a_ij = generate_polynomial(*Q,D);
                if Aij[[i,j]] == Rq::new(vec![]) {
                    Aij[[i,j]] = a_ij.clone();
                    Aij[[j,i]] = a_ij;
                }
            }
        }
        println!("Generated Aij!");
        /*
        for row in Aij.rows() {
            for poly in row {
                println!("{}", poly);
            }
        }
        */

        // random generation of random polynomial matrix Phi
        let mut Phi = Array2::from_elem((N,R), Rq::new(vec![])); 
        for i in 0..R {
            for j in 0..N {
                let phi_ji = generate_polynomial(*Q,D);
                Phi[[j,i]] = phi_ji;
            }
        }

        // Next, we need to compute 'b' such that this relation is equal to zero
        // So we compute the values of the relation and work backwards
        let mut a_product : Rq = Rq::new(vec![]);
        for i in 0..R {
            for j in 0..R {
                let vec = &S.column(i).to_vec();
                let vec2 = &S.column(j).to_vec();
                for k in 0..vec.len() {
                    //println!("VEC I ENTRY: {} ", vec[k]);
                    //println!("VEC J ENTRY: {} ", vec2[k]);
                    //println!("SANITY CHECK::: {} ", &S[[i,j]]);
                }

                let inner_prod = polynomial_vec_inner_product(&S.column(i).to_vec(), &S.column(j).to_vec());
                let prod = Aij[[i,j]].clone() * inner_prod;
                a_product = a_product.clone() + prod;
            }
        }

        let mut phi_product : Rq = Rq::new(vec![]);
        for i in 0..R {
            let inner_prod = polynomial_vec_inner_product(&Phi.column(i).to_vec(), &S.column(i).to_vec());
            phi_product = phi_product.clone() + inner_prod;
        }

        let b : Rq = &a_product + &phi_product;
        //println!("Generated b!");
        //println!("{}\n\n", b);

        //println!("A product: {}\n", a_product);
        //println!("Phi product: {}", phi_product);

        (Phi, Aij, b)
    }




    pub fn new(S: &Array2<Rq>) -> Self {
        let mut phi_k : Vec<Array2<Rq>> = vec![];
        let mut a_k : Vec<Array2<Rq>> = vec![];
        let mut b_k : Vec<Rq> = vec![];

        let mut phi_prime_k : Vec<Array2<Rq>> = vec![];
        let mut a_prime_k : Vec<Array2<Rq>> = vec![];
        let mut b_prime_k : Vec<Zq> = vec![];

        for k in 0..K {
            let phi : Array2<Rq>;
            let a : Array2<Rq>;
            let b : Rq;

            let (phi, a, b) = Self::gen_f(S);

            // TODO obviously cloning is not great here, but references are an even worse option
            // because no immediate way to do this without creating dangling pointer / not
            // compiling.. so we will potentially refactor later once minimally working.
            phi_k.push(phi.clone());
            a_k.push(a.clone());
            b_k.push(b.clone());

            phi_prime_k.push(phi);
            a_prime_k.push(a);
            b_prime_k.push(b.eval(Zq::from(0)) );
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



