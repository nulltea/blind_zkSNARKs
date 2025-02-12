# Verifiable Decryption for GBFV

We implemented the proof-of-decryption protocol presented in the following [paper](https://eprint.iacr.org/2024/1684) in C. We leveraged basic primitives used in [Lazer](https://eprint.iacr.org/2024/1846),
a library for lattice-based zero-knowledge proofs, and thoroughly extended it to construct our proof of decryption for [GBFV](https://eprint.iacr.org/2024/1587) ciphertexts.  


### Main source code
The main source code related to verifiable decryption can be found in:  

- Folder `vdec`
- Folder `scripts`
