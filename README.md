# Blind zkSNARKs
This repository contains all codes and assets regarding blind zkSNARKs. You can find:
1. **Proof of Decryption:** We implemented the proof-of-decryption protocol presented in our [paper] (https://eprint.iacr.org/2024/1684) in C. We leveraged basic primitives used in [Lazer](https://eprint.iacr.org/2024/1846), a library for lattice-based zero-knowledge proofs, and thoroughly extended it to construct our proof of decryption for [GBFV](https://eprint.iacr.org/2024/1587) ciphertexts.
   
3. **Performance estimation for Blind Fractal ([Magma](http://magma.maths.usyd.edu.au/magma/) scripts)**


### Proof of Decryption
The main source code related to verifiable decryption can be found in:  

- Folder `proof-of-decryption/vdec`
- Folder `proof-of-decryption/scripts`


### Team:
* Anonymous for now


### Acknowledgement
We thank Robin Geelen for helping us use the GBFV implementation and Vadim Lyubashevsky and Patrick Steuer for helping us use the Lazer library. A large part of the codebase is copied verbatim from [Lazer](https://eprint.iacr.org/2024/1846), produced by IBM Research Europe.
