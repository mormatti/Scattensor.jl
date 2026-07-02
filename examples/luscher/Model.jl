# =============================================================================
# Model.jl  --  Model definition and library-isolation wrappers
# =============================================================================
#
# Periodic transverse-field Ising model with a longitudinal field:
#
#     H = -J Σ_j Z_j Z_{j+1}  -  h_x Σ_j X_j  -  h_z Σ_j Z_j ,     Z_L ≡ Z_0
#
# with PBC, local dimension d = 2.  X = σ^x, Z = σ^z.
#
# IMPORTANT (design requirement from the assignment):
#   Every call into the underlying `Scattensor` library is isolated behind a
#   clearly named wrapper function in this file.  If you want to swap the
#   library (e.g. use a tensor-network backend, or your own ED routines),
#   you only need to re-point the three wrappers
#       build_hamiltonian, build_translation, raw_dispersion
#   and nothing else in the pipeline changes.
# =============================================================================

using SparseArrays
using LinearAlgebra

# ----------------------------------------------------------------------------
# Parameter struct (assignment-specified)
# ----------------------------------------------------------------------------
struct IsingParams
    J::Float64
    hx::Float64
    hz::Float64
end

IsingParams(; J = 1.0, hx = 0.6, hz = 0.05) = IsingParams(J, hx, hz)

# ----------------------------------------------------------------------------
# XXZ spin-1/2 chain (Bethe-ansatz integrable for ALL Δ):
#
#     H = -J Σ_j [ Sx_j Sx_{j+1} + Sy_j Sy_{j+1} + Δ Sz_j Sz_{j+1} ] - h Σ_j Sz_j
#
# Easy-axis ferromagnet (Δ > 1): the exact product state |↑↑…↑⟩ is the vacuum,
# single magnons have the exact dispersion ε(k)=J(Δ−cos k)+h, and the two-magnon
# S-matrix is known in closed form.  The uniform field h shifts an n-magnon
# state by n·h; choosing h > J(3−Δ) lifts the whole two-magnon sector above the
# one-magnon band (isolating a clean dispersion) WITHOUT changing the scattering
# phase (the field commutes with H and acts trivially within each Sz sector).
# ----------------------------------------------------------------------------
struct XXZParams
    J::Float64
    Δ::Float64
    h::Float64
end

XXZParams(; J = 1.0, Δ = 1.5, h = 2.0) = XXZParams(J, Δ, h)

# ----------------------------------------------------------------------------
# 1D Hubbard model (Lieb–Wu integrable), local dimension d = 4:
#
#     H = -t Σ_{j,σ}(c†_{jσ}c_{j+1,σ}+h.c.) + U Σ_j n_{j↑}n_{j↓}
#         + μ Σ_j(n_{j↑}+n_{j↓}) + Λ S²_tot P_{N=2}
#
# We benchmark the TWO-ELECTRON sector over the EMPTY vacuum:
#   * μ makes |0…0⟩ the ground state (one electron costs μ-2t>0);
#   * the (1↑,1↓) sector splits into a spin SINGLET (sees U) and a TRIPLET
#     (spatially antisymmetric ⇒ free).  The penalty Λ S²_tot P_{N=2} lifts the
#     triplet (S=1) away, isolating the interacting singlet — without touching
#     the vacuum or the one-electron band (P_{N=2} restricts it to N=2).
# Because each spin species has ≤1 particle in these sectors there are NO
# fermionic signs (Jordan–Wigner strings act trivially), so a hard-core-boson
# representation is exact here.
# ----------------------------------------------------------------------------
struct HubbardParams
    t::Float64
    U::Float64
    μ::Float64      # chemical potential: empty state = vacuum
    Λ::Float64      # S²_tot penalty: isolates the two-electron spin singlet
end

HubbardParams(; t = 1.0, U = 4.0, μ = 3.0, Λ = 10.0) = HubbardParams(t, U, μ, Λ)

# local Hilbert-space dimension of each model
local_dim(::IsingParams) = 2
local_dim(::XXZParams)   = 2
local_dim(::HubbardParams) = 4

# ----------------------------------------------------------------------------
# Run-configuration struct (assignment-specified)
# ----------------------------------------------------------------------------
struct LuscherConfig
    Ls::Vector{Int}                 # system sizes to diagonalize
    Ksectors::Vector{Int}           # total-momentum sectors Kint (K = 2π Kint / L)
    nlevels::Int                    # eigenstates per momentum sector
    rmax_dispersion::Int            # number of cosine harmonics in the ε(k) fit
    max_two_particle_levels::Int    # how many two-particle levels to try to match per K
end

LuscherConfig(; Ls = [8, 10, 12, 14],
                Ksectors = [0],
                nlevels = 12,
                rmax_dispersion = 4,
                max_two_particle_levels = 6) =
    LuscherConfig(Ls, Ksectors, nlevels, rmax_dispersion, max_two_particle_levels)

# ----------------------------------------------------------------------------
# Local (two-site) Hamiltonian term.
#
# `summation_local(H0, d, L; pbc=true)` sums a two-site operator H0 over every
# bond of the ring.  Single-site fields are split symmetrically across the two
# sites of a bond (factor 1/2) so that, after summing over all L bonds, each
# site is counted exactly once.
#
#     H0 = -J (Z⊗Z)  -  h_x/2 (X⊗I + I⊗X)  -  h_z/2 (Z⊗I + I⊗Z)
# ----------------------------------------------------------------------------
function ising_local_term(params::IsingParams)
    σx = SparseMatrixCSC([0.0 1.0; 1.0 0.0])
    σz = SparseMatrixCSC([1.0 0.0; 0.0 -1.0])
    id = SparseMatrixCSC([1.0 0.0; 0.0 1.0])
    ⊗(A, B) = kron(A, B)

    ZZ = σz ⊗ σz
    Xsym = 0.5 * (σx ⊗ id + id ⊗ σx)
    Zsym = 0.5 * (σz ⊗ id + id ⊗ σz)

    return -params.J * ZZ - params.hx * Xsym - params.hz * Zsym
end

# Two-site local term for the XXZ chain (complex Hermitian; Sα = σα/2).
function xxz_local_term(params::XXZParams)
    σx = SparseMatrixCSC(ComplexF64[0 1; 1 0])
    σy = SparseMatrixCSC(ComplexF64[0 -im; im 0])
    σz = SparseMatrixCSC(ComplexF64[1 0; 0 -1])
    id = SparseMatrixCSC(ComplexF64[1 0; 0 1])
    ⊗(A, B) = kron(A, B)

    XX = (σx ⊗ σx) / 4
    YY = (σy ⊗ σy) / 4
    ZZ = (σz ⊗ σz) / 4
    Zsym = 0.5 * (σz ⊗ id + id ⊗ σz) / 2        # (1/2)(Sz⊗I + I⊗Sz)

    return -params.J * (XX + YY + params.Δ * ZZ) - params.h * Zsym
end

# Hubbard d=4 site operators on |0>,|↑>,|↓>,|↑↓>  (hard-core boson rep, no JW).
function _hubbard_site_ops()
    id = SparseMatrixCSC{ComplexF64}(I, 4, 4)
    nup = spdiagm(0 => ComplexF64[0, 1, 0, 1])
    ndn = spdiagm(0 => ComplexF64[0, 0, 1, 1])
    bup = sparse(ComplexF64[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0])  # ↑: |↑>→|0>, |↑↓>→|↓>
    bdn = sparse(ComplexF64[0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0])  # ↓: |↓>→|0>, |↑↓>→|↑>
    Sz = 0.5 * spdiagm(0 => ComplexF64[0, 1, -1, 0])
    Sp = sparse(ComplexF64[0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0])    # S+ = |↑><↓|
    Sm = sparse(ComplexF64(1) * Sp')
    return (; id, nup, ndn, bup, bdn, Sz, Sp, Sm)
end

# Two-site term: hopping (both spins) + on-site U and μ split symmetrically.
function hubbard_local_term(p::HubbardParams)
    o = _hubbard_site_ops(); ⊗(A, B) = kron(A, B)
    sl(op) = 0.5 * (op ⊗ o.id + o.id ⊗ op)
    hop = -p.t * (o.bup' ⊗ o.bup + o.bup ⊗ o.bup' + o.bdn' ⊗ o.bdn + o.bdn ⊗ o.bdn')
    return hop + sl(p.U * (o.nup * o.ndn)) + sl(p.μ * (o.nup + o.ndn))
end

# Isolation operators for the two-electron benchmark (all DIAGONAL ⇒ cheap, and
# commuting with the translation T):
#   * (Sᶻ_tot)² P_{N=2}  lifts the 2↑ / 2↓ states (total Sᶻ=±1) within N=2,
#                        keeping the 1↑1↓ sector (Sᶻ=0);
#   * P_{N≥3}            lifts all ≥3-particle states (which otherwise overlap
#                        the two-particle energy window).
# The remaining 1↑1↓ sector still contains both the spin singlet (spatially
# symmetric, sees U) and triplet (antisymmetric, free); they are separated by
# REFLECTION parity at K=0 (handled in `compute_momentum_spectrum`).  NB: the
# hard-core-boson representation does not carry the Hubbard SU(2), so S²_tot is
# NOT a usable quantum number here — reflection parity is the right separator.
function hubbard_isolation_ops(L::Int)
    o = _hubbard_site_ops(); ⊗(A, B) = kron(A, B)
    sl(op) = 0.5 * (op ⊗ o.id + o.id ⊗ op)
    Sz = summation_local(sl(o.Sz), 4, L; pbc = true)
    N  = summation_local(sl(o.nup + o.ndn), 4, L; pbc = true)
    diagN = real.(diag(N)); diagSz = real.(diag(Sz))
    PN3 = spdiagm(0 => ComplexF64[diagN[i] > 2.5 ? 1.0 : 0.0 for i in eachindex(diagN)])
    Sz2PN2 = spdiagm(0 => ComplexF64[(abs(diagN[i] - 2) < 1e-6 ? diagSz[i]^2 : 0.0) for i in eachindex(diagN)])
    return Sz2PN2, PN3
end

# Generic accessor — add a method here to support a new model.
local_hamiltonian_term(p::IsingParams) = ising_local_term(p)
local_hamiltonian_term(p::XXZParams)   = xxz_local_term(p)
local_hamiltonian_term(p::HubbardParams) = hubbard_local_term(p)

# ============================================================================
#  LIBRARY-ISOLATION WRAPPERS  (the only place Scattensor is touched)
# ============================================================================

"""
    build_hamiltonian(L, params) -> AbstractMatrix

Return the finite periodic-boundary-condition Hamiltonian on `L` sites.
Works for any model that defines `local_hamiltonian_term(params)`.
Wrapper around `Scattensor.summation_local`.
"""
function build_hamiltonian(L::Int, params)
    H0 = local_hamiltonian_term(params)
    return summation_local(H0, local_dim(params), L; pbc = true)   # <-- library call
end

# Hubbard adds the (2↑/2↓ + ≥3-particle) isolation penalties.
function build_hamiltonian(L::Int, params::HubbardParams)
    H = summation_local(hubbard_local_term(params), 4, L; pbc = true)
    Sz2PN2, PN3 = hubbard_isolation_ops(L)
    return H + params.Λ * Sz2PN2 + 100.0 * PN3     # 100 ≫ two-particle bandwidth
end

"""
    build_translation(L, params) -> AbstractMatrix

Return the one-site lattice translation operator T (unitary, [H,T]=0) for the
model's local dimension.  Wrapper around `Scattensor.operator_translation`.
"""
function build_translation(L::Int, params)
    return operator_translation(SparseMatrixCSC, local_dim(params), L)  # <-- library call
end

"""
    reflection_operator(L, params) -> Union{AbstractMatrix, Nothing}

Spatial-reflection operator used to separate scattering channels at K=0.
Returned only for models that need it (Hubbard: even = interacting singlet,
odd = free triplet); `nothing` otherwise.
"""
reflection_operator(L::Int, params) = nothing
reflection_operator(L::Int, params::HubbardParams) =
    operator_reflection(SparseMatrixCSC, 4, L)   # <-- library call

"""
    raw_dispersion(H, T, L; nlevels) -> Vector{BlochState}

Diagonalize `H` in each momentum sector and return Bloch states carrying
`energy` and `koverpi` (= k/π).  Wrapper around `Scattensor.dispersion_relation`.

This is the single point through which momentum-resolved spectra enter the
pipeline.  Swap it to change the diagonalization backend.
"""
function raw_dispersion(H, T, L::Int; nlevels::Int)
    return dispersion_relation(H, T, L; nlevels = nlevels)   # <-- library call
end
