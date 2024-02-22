use proptest::prelude::*;
use proptest::test_runner::{Config, TestRunner};
use labrador_snark::constants::*;
use labrador_snark::util::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1))]

    #[test]
    fn test_JL_norm_preservation(v in prop::collection::vec(0..256i128, N*(D as usize))) {
        let original_norm = l2_norm(&v);
        let jl_norm = l2_norm(&jl_project_gen(&v));
        prop_assert!(jl_norm <= 11.3137*original_norm, "Projection norm: {}. Original L2 norm: {}", jl_norm, original_norm);
    }
}



