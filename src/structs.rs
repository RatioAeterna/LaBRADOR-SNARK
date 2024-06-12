use ndarray::Array2;
use std::collections::HashMap;
use std::mem;
use rand::Rng;
use rand::RngCore;
use rand_chacha::ChaCha20Rng;
use rand_core::SeedableRng as RandCoreSeedableRng;
use num_bigint::BigUint;

use crate::algebraic::*;
use crate::constants::*;
use crate::util::*;


// NOTE: Consider that for now, whoever is sampling from the CRS has to sample everything in a
// strict ordering..
// I.e., you sample the elements of "A" once, then the elements of "B", etc. every polynomial
// sample is going to increase the offset from the base seed in a predictable way until you reset it.
pub struct CRS<'a> {
    // Difference between these two: 'offset_seed' is the running counter that we keep
    base_seed: [u8; 32],
    constants: &'a RuntimeConstants
}

impl<'a> CRS<'a> {

    pub fn random_oracle_gen(&self, offset_seed : &mut [u8; 32]) -> Rq {
        let mut coefficients : Vec<Zq> = Vec::new();

        for _ in 0..D {
            coefficients.push(self.generate_random_coeff(offset_seed));
            //offset_seed = self.increment_seed(&offset_seed);
            self.increment_seed(offset_seed);
        }
        // convert coefficient vector to Rq and return
        Rq::new(coefficients)
    }

    fn biguint_to_bytes(bigint: BigUint) -> [u8; 32] {
        let mut bytes = [0u8; 32];
        let bigint_bytes = bigint.to_bytes_be();
        let start_index = 32 - bigint_bytes.len();
        bytes[start_index..].copy_from_slice(&bigint_bytes);
        bytes
    }

    pub fn fetch_A_row(&self, row : usize) -> Vec<Rq> {
        // first, set the seed properly

        // Convert base_seed to BigUint using big-endian byte order
        let base_seed_bigint = BigUint::from_bytes_be(&self.base_seed);

        // Calculate the increment
        let increment = BigUint::from(row) * BigUint::from(self.constants.N * D);

        // Add the increment to the base seed
        let offset_seed_bigint = base_seed_bigint + increment;

        // Convert back to big-endian byte array and create offset_seed
        let mut offset_seed = Self::biguint_to_bytes(offset_seed_bigint);

        // lastly, fetch and return the row
        self.fetch_next_n(self.constants.N, &mut offset_seed)
    }

    pub fn fetch_B_ik_row(&self, i : usize, k : usize, row : usize) -> Vec<Rq> {
        let base_seed_bigint = BigUint::from_bytes_be(&self.base_seed);

        // add offset from A
        let offset_A = BigUint::from(self.constants.KAPPA * self.constants.N * D);

        // get the size of a SINGLE B matrix
        let size_B_mat = BigUint::from(self.constants.KAPPA_1 * self.constants.KAPPA);
        let increment = offset_A + BigUint::from(i * self.constants.T_1 as usize + k)*size_B_mat + BigUint::from(row)*BigUint::from(self.constants.KAPPA * D);

        let offset_seed_bigint = base_seed_bigint + increment;
        let mut offset_seed = Self::biguint_to_bytes(offset_seed_bigint);

        self.fetch_next_n(self.constants.KAPPA, &mut offset_seed)
    }

    pub fn fetch_C_ijk(&self, i : usize, j : usize, k : usize) -> Vec<Rq> {
        let base_seed_bigint = BigUint::from_bytes_be(&self.base_seed);

        // add offset from A
        let offset_A = BigUint::from(self.constants.KAPPA * self.constants.N * D);

        // also add offset from all B matrices
        let size_B_mat = BigUint::from(self.constants.KAPPA_1 * self.constants.KAPPA);
        let num_B_matrices = BigUint::from(self.constants.R * self.constants.T_1 as usize);

        // figure out how far into C we are..
        let offset_C = BigUint::from(k + (self.constants.T_1 as usize)*(i*self.constants.R - i*(i-1)/2 + (j-i)));
        let increment = offset_A + num_B_matrices*size_B_mat*D + offset_C*BigUint::from(self.constants.KAPPA_2*D);

        let offset_seed_bigint = base_seed_bigint + increment;
        let mut offset_seed = Self::biguint_to_bytes(offset_seed_bigint);

        self.fetch_next_n(self.constants.KAPPA_2, &mut offset_seed)
    }

    pub fn fetch_D_ijk(&self, i : usize, j : usize, k : usize) -> Vec<Rq> {
        let base_seed_bigint = BigUint::from_bytes_be(&self.base_seed);

        // add offset from A
        let offset_A = BigUint::from(self.constants.KAPPA * self.constants.N * D);

        // also add offset from all B matrices
        let size_B_mat = BigUint::from(self.constants.KAPPA_1 * self.constants.KAPPA);
        let num_B_matrices = BigUint::from(self.constants.R * self.constants.T_1 as usize);

        // also add offset from all of C column vectors 
        let num_C_matrices = BigUint::from(self.constants.R * (self.constants.R + 1)/2);

        // figure out how far into D we are...
        let offset_D = BigUint::from(k + (self.constants.T_1 as usize)*(i*self.constants.R - i*(i-1)/2 + (j-i)));
        let increment = offset_A + num_B_matrices*size_B_mat*D + num_C_matrices*BigUint::from(self.constants.KAPPA_2*D) + offset_D*BigUint::from(self.constants.KAPPA_2*D);

        let offset_seed_bigint = base_seed_bigint + increment;
        let mut offset_seed = Self::biguint_to_bytes(offset_seed_bigint);

        self.fetch_next_n(self.constants.KAPPA_2, &mut offset_seed)
    }


    fn fetch_next_n(&self, n : usize, offset_seed : &mut [u8; 32]) -> Vec<Rq> {
        let mut res : Vec<Rq> = Vec::with_capacity(n);
        for _ in 0..n {
            res.push(self.random_oracle_gen(offset_seed));
        }
        res
    }

    fn increment_seed<'b>(&self, offset_seed: &'b mut [u8; 32]) -> &'b mut [u8; 32] {
        for byte in offset_seed.iter_mut().rev() {
            if *byte == 255 {
                *byte = 0;
            } else {
                *byte += 1;
                break;
            }
        }
        offset_seed
    }

    fn generate_random_coeff(&self, offset_seed : &[u8; 32]) -> Zq {
        //let offset_seed_bytestring = Self::biguint_to_bytes(*offset_seed);
        let mut rng = ChaCha20Rng::from_seed(*offset_seed); // Initialize RNG with a 256-bit seed
        Zq::from(rng.gen_range(0..*Q))
    }

    pub fn new(constants: &'a RuntimeConstants) -> Self {

        let mut rng = rand::thread_rng();
        let mut base_seed = [0u8; 32];
        rng.fill(&mut base_seed); // Generate a 256-bit base seed
                                  
        if is_verbose() {
            println!("Generated CRS random PRG seed... {:?}", base_seed);
        }

        let mut offset_seed = base_seed;

        CRS {
            base_seed,
            constants,
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

        if is_verbose() {
            println!("Generated phi!");
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
        if is_verbose() {
            println!("Generated b!");
        }

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
