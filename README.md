# Scattensor.jl

[![Build Status](https://github.com/mormatti/Scattensor.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mormatti/Scattensor.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mormatti/Scattensor.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mormatti/Scattensor.jl)

This package has the goal to provide functions and algorithms to simulate quasi-particle scattering in one-dimensional uniform discrete quatum many-body systems using Tensor Network methods (i.e. Matrix Product States).

### Assumptions (for the moment)

Let's consider a quantum many-body system with Hilbert space $\mathcal H$ and Hamiltonian $\hat H$. We require that

- The system 1-dimensional (space dimension);
- The system is discrete, on a (equally spaced) lattice.
- The system has a uniform local dimension $d$;
- For PBC of the system, it is defined a translation operator $\hat T$ (lattice geometry);
- $\hat H$ is translational invariant, i.e. $[\hat H, \hat T]=0$;
- The system is gapped (or, at least the gap is $\epsilon > 0$);
- The groundstate is non-degenerate (no SSB).

From these, some inputs:

- The local dimension of the system $d$ (of both ℋ and ℋ′);
- The number of sites L of the small PBC system ℋ′;
- The basis of the groundstate sector of ℋ′;
- The basis of the particle sector of ℋ′;
- The local Hamiltonian Hⱼ in matrix form (p-local operator like);
- The translation operator T for the two sectors;
- The translation operator $T$

We put the following limits:

- For $\dim \mathcal H < 20000$ we use exact diagonalization;
- For higher values of $\dim \mathcal H$  we use the Arnoldi method;
- For $d > 10$ we should use purely TN methods.

## Installation

Work in progress.

## Documentation

Work in progress.

## Examples

Work in progress.
