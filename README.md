# LaBRADOR-SNARK
Robust implementation of the LaBRADOR Post-Quantum Cryptographic Proof Protocol.

This is still a work-in-progress, though it's nearing completion.


## TODO
- [ ] Proper Generation of the integer modulus Q, selecting an appropriate prime so we're sure we're using the ring Rq
- [ ] Remove *all* clones, etc., basic optimization
- [ ] Fix floating point precision error with Check 14 that occurs with high values of Q
- [ ] Finish Recursion implementation, ideally benchmark with that 
- [ ] Fiat-Shamir implementation for non-interactivity
- [ ] Full ZK implementation
