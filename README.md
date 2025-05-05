# Scattensor.jl

[![Build Status](https://github.com/mormatti/Scattensor.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mormatti/Scattensor.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mormatti/Scattensor.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mormatti/Scattensor.jl)

This package has the goal to provide functions and algorithms to simulate quasi-particle scattering in one-dimensional uniform discrete quatum many-body systems using Tensor Network methods (i.e. Matrix Product States).

### Assumptions (for the moment)

Let's consider a quantum many-body system with Hilbert space $\mathcal H$ and Hamiltonian $H$. We require that

- **1D** **space** dimension: the space dimension of the system is $D=1$. Next step go higher dimensions $D>1$;
- **Discretization** of the space: The system is discrete, on lattice, with constant lattice spacing $a=1$. For some systems like Lattice Gauge Theories, the goal is to get near the continuum limit $a \to 0$.
- **Uniformity** of local dimension: The system has a uniform local dimension $d$. This also means that for PBC it is always defined a translation operator $\hat T$ of a unit site (on the right); Next step: generalize to the **Periodicity** $n$ of the local dimension (so that the translation operator is still defined). An example could be staggered fermions.
- **Translation invariance**: $H$ is translational invariant, i.e. $[H, T]=0$; Next steps involve generalizing this translational invariance to $T^n$ or a more generic translation.
- **Reflection** invariance: there exist a reflection operator (i.e. an operator such that $TR=RT^\dag$) under which the Hamiltonian is invariant (i.e. $[H,R] = 0$)
- **Locality** of the interaction: the Hamiltonian $H$ can be written as a sum of local Hamiltonians $H = \sum_{j}H_j$, where $H_j = T^j H_0 T^{j\dag}$. $H_0$ is called local Hamiltonian and by definition has a limited support which does not scale with the system size $L$.
- **Gapped** Hamiltonian: or, at least, the gap $\Delta > 0$  is not near the machine precision; The limitation to tackle $\Delta = 0$ is also due to Tensor Network area-law theorems. However, nothing forbids to try to simulate systems with $\Delta = 0$ in an approximate way.
- **Non-degeneration** of the groundstate (no SSB). Next step: include SSB phases (kink scattering).

## Installation

Work in progress.

## Documentation

Work in progress.

## Examples

Work in progress.
