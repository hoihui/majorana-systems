# majorana-systems
Simulation of non-interacting Majorana Hamiltonians

# Functions

* `H_1Ds.m`: generate Hamiltonian with general, non-uniform spin-orbin couplings, Zeeman splittings (in all three dimensions) and s-wave pairing
* `H_1Dp.m`: generate Hamiltonian with p-wave pairing.
* `H_1Ds_Fourier.m`: same as `H_1Ds.m` but in momentum space with `nband` parameter specifying number of bands to output. Used for trunction for large systems.
* `H_2Ds.m`: for 2D s-wave system.

More complete documentation can be found in each function's file.

# Examples
Pleaes execute section-by-section (separated by `##` in matlab) in SCp.m and SCs.m. Section headings explain the intentions of the code blocks.

* `SCs.m`: codes related to Hamiltonians with s-wave pairing.
* `SCp.m`: codes related to Hamiltonians with p-wave pairing.
