use proptest::prelude::*;
use proptest::test_runner::{Config, TestRunner};
use labrador_snark::constants::*;
use labrador_snark::algebraic::*;
use labrador_snark::util::*;
use std::sync::atomic::{AtomicBool, Ordering};

proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    #[test]
    fn test_JL_norm_preservation(v in prop::collection::vec(0..2^16i128, N*(D as usize))) {
        let original_norm = l2_norm(&v);
        let jl_norm = l2_norm(&jl_project_gen(&v));
        prop_assert!(jl_norm <= 2.0*11.3137*original_norm, "Projection norm: {}. Original L2 norm: {}", jl_norm, original_norm);
    }

    #[test]
    fn test_linearity_of_poly_vec_inner_product(a in prop::collection::vec(any::<Rq>(), N as usize),
                                                b in prop::collection::vec(any::<Rq>(), N as usize),
                                                c in any::<Zq>()) {
        NTT_ENABLED.store(false, Ordering::SeqCst);
        let product_ab = polynomial_vec_inner_product(&a, &b);
        println!("product ab! {}", product_ab);
        let product_ab_scaled = polynomial_vec_inner_product(&a, &scale_poly_vec_int(&b, &c));
        println!("scaled of ab by c! {}", scale_polynomial_int(&product_ab, &c));

        // Testing linearity: c * inner_product(a, b) should equal inner_product(a, c*b)
        prop_assert_eq!(product_ab_scaled, scale_polynomial_int(&product_ab, &c), "Linearity not preserved in polynomial vector inner products, for some reason...");
    }


    #[test]
    fn test_linearity_ntt_enabled(a in prop::collection::vec(any::<Rq>(), N as usize),
                                                b in prop::collection::vec(any::<Rq>(), N as usize),
                                                c in any::<Zq>()) {
        
        NTT_ENABLED.store(true, Ordering::SeqCst);
        let product_ab = polynomial_vec_inner_product(&a, &b);
        println!("product ab! {}", product_ab);
        let product_ab_scaled = polynomial_vec_inner_product(&a, &scale_poly_vec_int(&b, &c));
        println!("scaled of ab by c! {}", scale_polynomial_int(&product_ab, &c));

        // Testing linearity: c * inner_product(a, b) should equal inner_product(a, c*b)
        prop_assert_eq!(product_ab_scaled, scale_polynomial_int(&product_ab, &c), "Linearity not preserved in polynomial vector inner products, for some reason...");
    }


    #[test]
    fn sigma_inv_invariant(a in prop::collection::vec(any::<Zq>(), N*(D as usize)),
                           b in prop::collection::vec(any::<Zq>(), N*(D as usize))) {
        let inner_prod : Zq = vec_inner_product_Zq(&a,&b);
        let poly_vec_a : Vec<Rq> = concat_coeff_reduction(&a);
        let poly_vec_b : Vec<Rq> = concat_coeff_reduction(&b);
        let inv_a = sigma_inv_vec(&poly_vec_a);

        let poly_prod : Rq = polynomial_vec_inner_product(&inv_a, &poly_vec_b);
        //let poly_eval_to_int : i128 = i128::from(poly_prod.eval(Zq::from(0)));
        let poly_eval = poly_prod.eval(Zq::from(0));

        //prop_assert!(inner_prod == poly_eval_to_int, "Conjugation Autmorphism invariant failed! Regular inner product {} should equal constant term {}", inner_prod, poly_eval_to_int);
        prop_assert!(inner_prod == poly_eval, "Conjugation Autmorphism invariant failed! Regular inner product {} should equal constant term {}", inner_prod, poly_eval);
    }
}




/*
#[test]
fn util_test() {
    // test polynomial generation
    let bigval : Zq = Zq::new(4294967295);
    let one : Zq = Zq::one();
    let minusone : Zq = Zq::from(-1);
    println!("bigval: {}", bigval);
    println!("sum with one: {}", bigval + one);
    println!("-1 in Zq: {}", minusone);
    println!("product by two: {}", bigval*Zq::from(2));



    let mut data : Vec<Zq> = vec![Zq::zero(); 81];
    let mut data2 : Vec<Zq> = vec![Zq::zero(); 81];
    data[80] = Zq::from(4294967295 as i128);
    data2[1] = Zq::from(2);

    println!("Creating big poly 1!");
    let too_big_poly = Rq::new(data);
    println!("Creating big poly 2!");
    let too_big_poly2 = Rq::new(data2);
    println!("poly: {}", too_big_poly);
    println!("poly2: {}", too_big_poly2);

    let prod_poly = too_big_poly * too_big_poly2;
    println!("prod poly: {}", prod_poly);
    /*
    let p1 = generate_polynomial(Q, D);

    println!("Random polynomial:");
    println!("{}", p1);


    // test random sampling Zq
    let sample = random_sample_Zq(Q, 256);
    println!("String of integers mod q:");
    println!("{:?}", sample);

    // test scaling polynomial by a float < 1
    let scale: f32 = 0.5;
    let p2 = scale_polynomial(p1, scale);
    println!("Polynomial scaled by {:.32}:", scale);
    println!("{}", p2);

    let S = generate_witness();
    println!("Generated witness matrix S:");
    for row in S.rows() {
        for poly in row {
            println!("{}", poly);
        }
    }

    // test norm computation of Array2 of polynomials
    let witness_norm : f64 = compute_total_norm(S);
    println!("Computed norm of witness matrix: {:.64}", witness_norm);

    println!("Testing polynomial vec inner product:");
    let mut vec_a : Vec<Polynomial<i128>> = vec![];
    let mut vec_b : Vec<Polynomial<i128>> = vec![];
    for i in 0..5 {
        let mut new_poly = Polynomial::new(vec![0i128,1,2,3]);
        vec_a.push(new_poly.clone());
        vec_b.push(new_poly);
    }

    let product : Polynomial<i128> = Polynomial::new(vec![0i128,1,2,3]) * Polynomial::new(vec![0i128,1,2,3]) * Polynomial::new(vec![5i128]);

    let prod_poly = polynomial_vec_inner_product(vec_a, vec_b);
    assert!(prod_poly == product, "Polynomial inner products are broken.");
    println!("Polynomial products are working!");
    println!("{}", prod_poly);
    */
}
*/
