# Scattensor.jl

[![Build Status](https://github.com/mormatti/Scattensor.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mormatti/Scattensor.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mormatti/Scattensor.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mormatti/Scattensor.jl/branch/main)

**Scattensor.jl** is a Julia package for simulating quasi-particle scattering in one-dimensional uniform discrete quantum many-body systems using Tensor Network methods (Matrix Product States/MPOs).

## Overview

This package provides tools and algorithms to:
- Compute dispersion relations and band structures for translation-invariant Hamiltonians
- Extract ground states and excited states using tensor network methods
- Construct Wannier functions for maximally localized states
- Compute S-matrices for scattering processes
- Perform time evolution using TDVP (Time-Dependent Variational Principle)
- Work with both matrix and tensor network representations of quantum operators

The package bridges abstract mathematical operators with efficient ITensor implementations, supporting both dense/sparse matrix representations and tensor network (MPS/MPO) representations.

## Installation

**Note:** This package is currently in active development (version 1.0.0-DEV). Installation instructions are under construction.

For development installation, you can add the package locally:

```julia
using Pkg
Pkg.develop(path="/path/to/Scattensor")
```

## Quick Start

```julia
using Scattensor
using ITensors, ITensorMPS
using SparseArrays

# Define local Hamiltonian H0 and system parameters
d = 2  # local dimension
L = 10  # system length
H0 = ...  # your local Hamiltonian

# Build full Hamiltonian and translation operator
H = summation_local(H0, d, L; pbc = true)
T = operator_translation(SparseMatrixCSC, d, L)

# Compute dispersion relation
disprel = dispersion_relation(H, T, L; nlevels = 10)
plot_disprel(disprel)

# Extract ground state and first band
ω_gs = pop_groundstate!(disprel)
band = pop_firstband!(disprel)
```

## Features

### Implemented Features

- **Dispersion Relations**: Compute energy bands and momentum eigenstates for translation-invariant systems
- **Ground State Extraction**: Extract and manipulate ground states from dispersion relations
- **Band Structure Analysis**: Access and filter states by energy bands and momentum
- **Wannier Functions**: Compute maximally localized Wannier functions with symmetry constraints
- **Operator Construction**: Build translation, reflection, and identity operators for both matrix and tensor representations
- **Tensor Network Operations**: Convert between matrix and MPO/MPS representations
- **Time Evolution**: TDVP-based time evolution for MPS states
- **Local Expectation Values**: Compute expectation values of local operators
- **Visualization**: Plot dispersion relations and complex-valued data

### Under Construction

- **S-Matrix Computation**: S-matrix functions exist but require testing and debugging
- **Documentation**: Comprehensive API documentation is being developed
- **Test Coverage**: Test suite is minimal and needs expansion
- **Installation**: Package registration and installation instructions are pending

## System Requirements and Assumptions

The package is designed for quantum many-body systems satisfying the following assumptions:

### Current Assumptions

- **1D Space Dimension**: The system is one-dimensional ($D=1$). Higher dimensions are planned for future releases.
- **Discrete Lattice**: The system is discrete with constant lattice spacing $a=1$ (can approach continuum limit $a \to 0$ for lattice gauge theories).
- **Uniform Local Dimension**: All sites have the same local dimension $d$. Generalization to periodic local dimensions (e.g., staggered fermions) is planned.
- **Translation Invariance**: The Hamiltonian commutes with the translation operator, $[H, T] = 0$. Generalizations to $T^n$ or generic translations are planned.
- **Reflection Invariance**: A reflection operator $R$ exists such that $TR = RT^\dagger$ and $[H, R] = 0$.
- **Local Interactions**: The Hamiltonian can be written as $H = \sum_j H_j$ where $H_j = T^j H_0 T^{j\dagger}$ and $H_0$ has limited support.
- **Gapped Spectrum**: The Hamiltonian has a gap $\Delta > 0$ not near machine precision. Gapless systems can be approximated but may violate tensor network area-law assumptions.
- **Non-Degenerate Ground State**: The ground state is unique (no spontaneous symmetry breaking). SSB phases and kink scattering are planned for future releases.

### Planned Extensions

- Higher dimensional systems ($D > 1$)
- Periodic local dimensions (staggered fermions)
- Generalized translation invariance ($T^n$)
- Spontaneous symmetry breaking phases
- Kink scattering

## Examples

Example scripts are provided in the `examples/` directory:

- **Ising Model** (`examples/ising/`):
  - Dispersion relation computation
  - Wannier function construction
  - Ground state extraction
  - Time evolution
  - S-matrix computation (under construction)

- **Bose-Hubbard Model** (`examples/bosehubbard/`):
  - Dispersion relation computation
  - Wannier function construction
  - Time evolution

See the example files for detailed usage patterns.

## Architecture

The package is organized in a four-tier architecture:

1. **Abstract Functions** (`abstract_functions/`): Generic interfaces for operators independent of representation
2. **Extensions**: Dual implementations for matrices and tensors:
   - **Matrix Extensions** (`matrix_extensions/`): Dense/sparse matrix implementations
   - **ITensor Extensions** (`itensor_extensions/`): ITensor/MPO/MPS implementations
3. **Core Functions** (`scattensor_functions/`): High-level physics algorithms
4. **Types & Utils** (`types/`, `utils/`): Type system and helper functions

## Dependencies

- **ITensors.jl / ITensorMPS.jl**: Tensor network backend
- **LinearAlgebra / SparseArrays**: Matrix operations
- **KrylovKit.jl**: Eigensolvers
- **Optim.jl**: Optimization
- **Plots.jl / Colors.jl**: Visualization
- **LaTeXStrings.jl**: LaTeX string support

See `Project.toml` for complete dependency list and version constraints.

## Contributing

Contributions are welcome! Please note that this package is in active development. When contributing:

- Follow the existing code structure and patterns
- Add tests for new functionality
- Update documentation as needed
- Note any limitations or TODOs in code comments

## License

This package is licensed under the MIT License. See `LICENSE` for details.

## Citation

If you use Scattensor.jl in your research, please cite:

```bibtex
@software{scattensor2023,
  author = {Morgavi, Mattia},
  title = {Scattensor.jl: Quasi-particle Scattering in 1D Quantum Many-Body Systems},
  year = {2023},
  url = {https://github.com/mormatti/Scattensor.jl}
}
```

## Status

**Current Version**: 1.0.0-DEV

This package is under active development. Core functionality for dispersion relations, ground states, and Wannier functions is implemented and functional. S-matrix computation and comprehensive documentation are in progress.

## Contact

For questions, issues, or contributions, please open an issue on [GitHub](https://github.com/mormatti/Scattensor.jl).
