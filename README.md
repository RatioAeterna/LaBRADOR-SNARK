# LaBRADOR-SNARK
Robust implementation of the LaBRADOR Post-Quantum Cryptographic Proof Protocol.

This is still a work-in-progress, though it's nearing completion.


## TODO
- [ ] Proper Generation of the integer modulus Q, selecting an appropriate prime so we're sure we're using the ring Rq
- [ ] Nit: computational verification of the algorithm used for the above point
- [ ] Investigate general slowness of implementation (perhaps expected)
- [ ] Remove *all* clones, etc., basic optimization
- [ ] Fix floating point precision error with Check 14 that occurs with high values of Q
- [ ] Finish Recursion implementation, ideally benchmark with that 
- [ ] Use lattice estimation tool for SIS to determine appropriate values for Q, D, KAPPA such that we have a security level of >= 2^128
- [ ] Fiat-Shamir implementation for non-interactivity
- [ ] Full ZK implementation
