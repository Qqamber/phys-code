# phy-codes

## DMRG

A finite(infinite) dmrg code based on White's paper.

Refer to [Density matrix formulation for quantum renormalization groups](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.2863) and [Density-matrix algorithms for quantum renormalization groups](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.10345) for details.

And also a dmrg code based on mps.

In mathematica, we use `TensorContract`, `TensorProduct`, `ArrayReshape` and `Transpose` to handle the index implementations. And using `Inactive`, `Activate` with `TensorContract`, `TensorProduct` can improve the performance of the contraction of tensors.

Refer to [Empter's note](https://github.com/empter/DMRGwithMPSandMPO/blob/master/DMRG_MPS.pdf) for the index convention in the code.

## iTEBD

A mathematica code for itebd and an itensor version based on Vidal's paper.

Refer to [Efficient classical simulation of slightly entangled quantum computations](https://arxiv.org/abs/quant-ph/0301063), [Efficient simulation of one-dimensional quantum many-body systems](https://arxiv.org/abs/quant-ph/0310089) and [Classical simulation of infinite-size quantum lattice systems in one spatial dimension](https://arxiv.org/abs/cond-mat/0605597) for details.

The itensor version code is based on [Itensor library v3](http://itensor.org/). You can refer to itensor's home page for the installation. A good tutorial on windows OS can refer to [Empter's note](https://github.com/empter/mps-tutorial), too.

## Model Check

A code aims to check the numerical results of the algorithms above for small spin systems, which has an itensor alike interface to handle the Hamiltonian terms.
