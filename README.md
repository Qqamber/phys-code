# phy-code-open

## DMRG

A finite and infinite dmrg based on White's paper.

Refer to [Density matrix formulation for quantum renormalization groups](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.2863) and [Density-matrix algorithms for quantum renormalization groups](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.10345) for details.

And also a dmrg code based on mps.

In mathematica, we use `TensorContract`, `TensorProduct`, `ArrayReshape` and `Transpose` to handle the index implementations.

Refer to [Empter's note](https://github.com/empter/DMRGwithMPSandMPO/blob/master/DMRG_MPS.pdf) for the index convention in the code.

## iTEBD

A mathematica code for itebd and an itensor version based on Vidal's paper.

Refer to [Efficient classical simulation of slightly entangled quantum computations](https://arxiv.org/abs/quant-ph/0301063), [Efficient simulation of one-dimensional quantum many-body systems](https://arxiv.org/abs/quant-ph/0310089) and [Classical simulation of infinite-size quantum lattice systems in one spatial dimension](Classical simulation of infinite-size quantum lattice systems in one spatial dimension) for details.

## Model Check

A code aims to check the numerical results of the algorithms above, which has an itensor alike interface to handle the Hamiltonian terms.