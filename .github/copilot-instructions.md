# Scattensor.jl Copilot Instructions

## Project Overview
**Scattensor.jl** is a Julia package for simulating quasi-particle scattering in 1D quantum many-body systems using Tensor Network Methods (Matrix Product States/MPOs). It bridges abstract mathematical operators with ITensor implementations for efficient quantum simulations.

## Core Architecture

### Module Organization (4-tier layering)
1. **Abstract Functions** (`abstract_functions/`): Generic interfaces for translation/reflection operators and Kronecker products independent of representation
2. **Extensions**: Dual implementations of abstract operations:
   - **Matrix Extensions** (`matrix_extensions/`): Dense/sparse matrix implementations
   - **ITensor Extensions** (`itensor_extensions/`): ITensor/MPO/MPS implementations
3. **Core Scattensor Functions** (`scattensor_functions/`): High-level physics algorithms (ground state, band structure, Wannier functions, S-matrices)
4. **Types & Utils** (`types/`, `utils/`): Type system and helper functions

### Key Design Patterns

#### Multi-Dispatch Type Hierarchy
The codebase heavily uses Julia's multiple dispatch with a generic/specific implementation pattern:
- **Generic interface** in `abstract_functions/operator_translation.jl`: `function operator_translation end`
- **Concrete implementations** split by representation (Matrix vs ITensor) with full method signatures
- Example: `operator_translation(MatrixType, d::Int, L::Int)` vs ITensor equivalents
- See `abstract_functions/` for the interface pattern; check both `matrix_extensions/` and `itensor_extensions/` for implementations

#### Representation-Agnostic Operators
Three types of operators are implemented twice (matrix & ITensor versions):
- **Translation** (`T`): Cyclic right-shift; satisfies translation invariance `[H, T] = 0`
- **Reflection** (`R`): Parity operator; satisfies `TR = RT†`
- **Identity**: Local identity operators for tensor construction

#### The `State{HilbertSpace, DataType}` Type
Custom generic type wrapping physics data with metadata:
```julia
mutable struct State{HilbertSpaceType, DataType}
    hilbspace::HilbertSpaceType
    data::DataType  # MPS/MPO or vector
    time::Real
end
```
- Enables polymorphic physics operations across different data representations
- `BlochState` is specialized for momentum eigenstates with energy/momentum fields

### Critical Functions (Entry Points)

| Function | File | Purpose | Key Inputs |
|----------|------|---------|-----------|
| `get_groundstate` | `scattensor_functions/get_groundstate.jl` | Extract lowest-energy state from vector | `Vector{<:BlochState}` |
| `mpo_from_matrix` | `itensor_extensions/mpo_from_matrix.jl` | Convert dense matrix → MPO for efficient tensors | `Matrix`, local dimension `d` |
| `wannier_symmetric` | `scattensor_functions/wannier_symmetric.jl` | Compute maximally-localized Wannier functions | Band states, parity op, inner product function |
| `apply_translation`/`apply_reflection` | `itensor_extensions/apply_*.jl` | Apply symmetry ops to MPS/MPO | State, operator |
| `partial_trace` | `itensor_extensions/partial_trace.jl` | Trace out subsystems (density matrix reduction) | MPS, sites to trace |

## Workflow Patterns

### Physics Simulation Pipeline
1. **Define Hamiltonian**: Build local Hamiltonian `H0` as matrix or MPO
2. **Compute Ground State**: Use DMRG (via ITensorMPS integration)
3. **Build Band Structure**: Compute Bloch states across k-points
4. **Analyze Operators**: Extract Wannier functions, S-matrices, expectation values
5. **Time Evolution**: Apply TDVP for dynamics

### When Adding Operator Functions
1. **Implement abstract interface first** in `abstract_functions/operator_*.jl`
2. **Add matrix version** in `matrix_extensions/operator_*.jl` (sparse or dense)
3. **Add ITensor version** in `itensor_extensions/operator_*.jl`
4. **Export from main module** in `src/Scattensor.jl`
5. See `operator_translation.jl` across all three folders as the canonical example

## External Dependencies & Integration

### Core Dependencies
- **ITensors.jl / ITensorMPS.jl**: Tensor network backend; provides `MPS`, `MPO`, `Index`, `ITensor`, `dmrg()`, `tdvp()`
- **LinearAlgebra / SparseArrays**: Dense/sparse matrix ops; used extensively in matrix_extensions
- **KrylovKit.jl**: Eigensolvers for band structure computation
- **Optim.jl**: Optimization (likely for Wannier function localization)
- **Plots.jl / Colors.jl**: Visualization of dispersion relations and complex-valued data

### Custom DMRG Observer
`custom_dmrg_observer.jl` hooks into ITensorMPS's DMRG loop for custom convergence tracking and observables collection during ground state search.

## Testing & Development

### Running Tests
```bash
julia> include("test/runtests.jl")
```
- Test infrastructure in place but minimal coverage (see `test/runtests.jl` — mostly precompilation validation)
- Examples in `examples/ising/` and `examples/bosehubbard/` serve as primary integration tests
- **Note**: Use `Revise.jl` for interactive development (loaded in examples and test file)

### Workflow Automation
- No build/CI pipeline visible in source (check GitHub Actions)
- Default tolerances set globally: `default_cutoff = 1e-10`, `default_maxdim = 100`
- Consider these when tuning numerical precision in new features

## Project Constraints & Assumptions

The Scattensor framework is built on restrictive assumptions (from README):
- **1D** lattice with constant spacing
- **Uniform local dimension** `d` across all sites
- **Translation invariant** Hamiltonian (generalizations future work)
- **Gapped** spectrum (area-law theorems depend on this)
- **Non-degenerate** ground state (no spontaneous symmetry breaking)

When extending the code, respect these assumptions or clearly document violations.

## Conventions

- **Complex numbers as eigenvector data**: Functions like `rgb_from_complex.jl` handle phase visualization
- **Rational momentum values**: Bloch states use `koverpi::Rational` (momentum as fraction of π)
- **Terminal output**: `cancel_terminal_line.jl` provides progress bar cleanup for long-running DMRG
- **Matrix wrapper constructors**: E.g., `SparseMatrixCSC([...])` for compatibility with sparse matrix API
- **Kron operator**: Custom `⊗` alias for `kron()` used in examples for readability

## Common Pitfalls

1. **Forgetting to export** new functions in `src/Scattensor.jl` — add to `include()` list with proper order (abstract → extensions → scattensor_functions)
2. **Mixing representations**: MPO functions won't work on dense matrices; ensure type consistency or convert via `mpo_from_matrix()`
3. **Mismatched local dimensions**: `d` must be consistent across all sites; validation in `is_uniform_localdim.jl` and similar
4. **Parity vs translation operator order**: Relations like `TR = RT†` are exact; check `apply_reflection.jl` for correct composition
