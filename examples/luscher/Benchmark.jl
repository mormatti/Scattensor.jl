# =============================================================================
# Benchmark.jl  --  exact / literature reference results
# =============================================================================
#
# Two independent analytic benchmarks for validating the Lüscher pipeline:
#
#  (A) FREE-FERMION point  (h_z = 0).  The transverse-field Ising chain is
#      Jordan-Wigner free fermions, with the EXACT single-particle dispersion
#
#          ε(k) = 2 √( J² + h_x² − 2 J h_x cos k ).
#
#      The bulk two-particle S-matrix is the constant S = −1 (pure statistics),
#      so the extracted phase shift should be flat (δ = const, mod π).  This
#      validates dispersion, momentum quantization and level bookkeeping.
#
#  (B) E8 point  (h_x = J critical, small h_z).  Ising CFT + magnetic field is
#      Zamolodchikov's integrable E8 theory.  Two clean references:
#        * mass ratios  m_n/m_1  (e.g. m_2/m_1 = golden ratio 1.6180…),
#        * the EXACT lightest-lightest elastic S-matrix
#              S_11(θ) = {2/3}{2/5}{1/15},
#          giving the elastic phase shift δ_11(θ).
#      Refs: Zamolodchikov, Int.J.Mod.Phys.A 4, 4235 (1989);
#            Dorey, hep-th/9810026; Coldea et al., Science 327, 177 (2010).
# =============================================================================

# ----------------------------------------------------------------------------
# (A) Exact free-fermion dispersion (h_z = 0)
# ----------------------------------------------------------------------------
"""
    exact_free_fermion_dispersion(params, k) -> Float64

Exact one-particle dispersion of the pure transverse-field Ising chain
ε(k) = 2√(J² + h_x² − 2 J h_x cos k).  Only meaningful at h_z = 0.
"""
function exact_free_fermion_dispersion(params::IsingParams, k::Real)
    return 2 * sqrt(params.J^2 + params.hx^2 - 2 * params.J * params.hx * cos(k))
end

# ----------------------------------------------------------------------------
# (B) E8 mass spectrum (ratios to the lightest mass m1)
# ----------------------------------------------------------------------------
const E8_MASS_RATIOS = [
    1.0,                # m1
    1.6180339887498949, # m2/m1 = 2cos(π/5)  (golden ratio φ)
    1.9890437907365871, # m3/m1 = 2cos(π/30)
    2.4048671724351026, # m4/m1
    2.9562952014788549, # m5/m1
    3.2183404585766740, # m6/m1
    3.8911568233405156, # m7/m1
    4.7833861168556960, # m8/m1
]

# ----------------------------------------------------------------------------
# (B) Exact E8 lightest-lightest S-matrix  S_11(θ) = {2/3}{2/5}{1/15}
#
# Block:  (sinh θ + i sin(π a)) / (sinh θ − i sin(π a)),  modulus 1 for real θ.
# Its phase is 2·atan( sin(π a) / sinh θ ).
# ----------------------------------------------------------------------------
const _E8_BLOCKS = (2 // 3, 2 // 5, 1 // 15)

"""
    e8_S11(θ) -> ComplexF64

Exact elastic S-matrix of two lightest E8 particles at relative rapidity θ.
"""
function e8_S11(θ::Real)
    s = sinh(θ)
    z = one(ComplexF64)
    for a in _E8_BLOCKS
        sa = sin(π * float(a))
        z *= (s + im * sa) / (s - im * sa)
    end
    return z
end

"""
    e8_delta11(θ) -> Float64

Elastic phase shift δ_11(θ) in the convention S = e^{2iδ}, i.e.
δ_11(θ) = Σ_a atan( sin(π a) / sinh θ ).  Goes from 3π/2 at θ→0⁺ to 0 at θ→∞.
"""
function e8_delta11(θ::Real)
    s = sinh(θ)
    return sum(atan(sin(π * float(a)) / s) for a in _E8_BLOCKS)
end

# ----------------------------------------------------------------------------
# Lattice ↔ rapidity map (K = 0 centre-of-mass frame)
#
# A single particle of lattice momentum q has energy ε(q) = m1 cosh β, so its
# rapidity is β = acosh(ε(q)/m1).  In the two-particle CM frame the relative
# rapidity is θ = 2β.  This uses the MEASURED dispersion ε(q) and the measured
# gap m1 = min ε, so it needs no separate velocity input.
# ----------------------------------------------------------------------------
"""
    rapidity_from_q(fit, q, m1) -> Float64

Relative rapidity θ = 2 acosh(ε(q)/m1) for a K=0 two-particle state with
relative lattice momentum q.
"""
function rapidity_from_q(fit::FourierDispersion, q::Real, m1::Real)
    r = epsilon(fit, q) / m1
    r = max(r, 1.0)                      # guard against tiny numerical r<1
    return 2 * acosh(r)
end

"""
    e8_phase_curve(fit, qs, m1) -> (θs, δs)

Map a grid of relative lattice momenta `qs` to E8 rapidities and return the
exact E8 elastic phase shift δ_11 on that grid (K=0).
"""
function e8_phase_curve(fit::FourierDispersion, qs, m1::Real)
    θs = [rapidity_from_q(fit, q, m1) for q in qs]
    δs = [e8_delta11(θ) for θ in θs]
    return θs, δs
end

# ============================================================================
#  (C) XXZ easy-axis ferromagnet — EXACT two-magnon Bethe-ansatz benchmark
#
#  Single magnon:   ε(k) = J (Δ − cos k) + h            (h = uniform field)
#  Two-magnon S-matrix at total momentum K=0 (relative momentum q), from the
#  Bethe ansatz (Bethe 1931; XXZ form e.g. Karbach–Müller, cond-mat/9809162):
#
#      S(q) = − (1 − Δ e^{−iq}) / (1 − Δ e^{iq}) ,     S = e^{2iδ}
#      δ(q) = π/2 + atan2( Δ sin q , 1 − Δ cos q ).
#
#  The field h shifts ε but leaves S and δ unchanged (non-interacting).  At
#  Δ = 0 (the XX point) this reduces to free fermions, S = −1, δ = π/2.
# ============================================================================

"""
    xxz_magnon_dispersion(params, k) -> Float64

Exact single-magnon dispersion ε(k) = J(Δ − cos k) + h of the XXZ ferromagnet.
"""
xxz_magnon_dispersion(params::XXZParams, k::Real) =
    params.J * (params.Δ - cos(k)) + params.h

"""
    xxz_two_magnon_S(params, K, q) -> ComplexF64
    xxz_two_magnon_S(params, q)    -> ComplexF64   (K = 0)

Exact two-magnon elastic S-matrix at total momentum K and relative momentum q.
From the coordinate Bethe ansatz, S(k₁,k₂)=−[e^{iK}−2Δe^{ik₁}+1]/[e^{iK}−2Δe^{ik₂}+1]
with k₁,₂=K/2±q, which factorizes to a function of cos(K/2):

    S_K(q) = −(cos(K/2) − Δ e^{−iq}) / (cos(K/2) − Δ e^{iq}).

At K=0 this is the headline formula −(1−Δe^{−iq})/(1−Δe^{iq}).
"""
xxz_two_magnon_S(params::XXZParams, K::Real, q::Real) =
    -(cos(K/2) - params.Δ * exp(-im * q)) / (cos(K/2) - params.Δ * exp(im * q))
xxz_two_magnon_S(params::XXZParams, q::Real) = xxz_two_magnon_S(params, 0.0, q)

"""
    xxz_two_magnon_delta(params, K, q) -> Float64
    xxz_two_magnon_delta(params, q)    -> Float64   (K = 0)

Exact two-magnon phase shift δ_K(q) = π/2 + atan2(Δ sin q, cos(K/2) − Δ cos q).
At K=0 it runs from 3π/2 at threshold to π/2 at q=π; raising K (lowering cos(K/2))
weakens the effective coupling, flattening δ toward π/2.
"""
xxz_two_magnon_delta(params::XXZParams, K::Real, q::Real) =
    π/2 + atan(params.Δ * sin(q), cos(K/2) - params.Δ * cos(q))
xxz_two_magnon_delta(params::XXZParams, q::Real) = xxz_two_magnon_delta(params, 0.0, q)

# ============================================================================
#  (D) 1D Hubbard model — EXACT two-electron (Lieb–Wu) benchmark
#
#  Empty vacuum; one electron is free: ε(k) = −2t cos k (+ chemical potential μ).
#  Two electrons of opposite spin form a spin SINGLET (symmetric in space → sees
#  U) and a spin TRIPLET (antisymmetric → free, δ=0).  The interacting singlet
#  two-body S-matrix at total momentum K=0 (relative momentum q) is, from the
#  Bethe ansatz / direct lattice two-body solution,
#
#      δ(q) = π/2 + atan( 4 t sin q / U ),   S(q) = −(U + 4 i t sin q)/(U − 4 i t sin q).
#
#  Limits: U→0 ⇒ δ=π (S=+1, free); U→∞ ⇒ δ=π/2 (S=−1, fermionized hard core).
#  Ref: Lieb & Wu, PRL 20, 1445 (1968); Essler et al., "The 1D Hubbard Model".
# ============================================================================

"""
    hubbard_dispersion(params, k) -> Float64

One-electron dispersion ε(k) = −2t cos k + μ.
"""
hubbard_dispersion(params::HubbardParams, k::Real) =
    -2 * params.t * cos(k) + params.μ

"""
    hubbard_singlet_S(params, K, q) -> ComplexF64
    hubbard_singlet_S(params, q)    -> ComplexF64   (K = 0)

Exact opposite-spin (singlet) two-electron S-matrix.  The interaction sees the
combination sin k₁ − sin k₂ = 2 cos(K/2) sin q, so

    S_K(q) = −(U + 4 t cos(K/2) sin q · i)/(U − 4 t cos(K/2) sin q · i).

(NB: this benchmark is only *used* at K=0 — the singlet/triplet separation relies
on K=0 reflection parity, which is not a good quantum number at K≠0.)
"""
hubbard_singlet_S(params::HubbardParams, K::Real, q::Real) =
    -(params.U + 4im * params.t * cos(K/2) * sin(q)) /
     (params.U - 4im * params.t * cos(K/2) * sin(q))
hubbard_singlet_S(params::HubbardParams, q::Real) = hubbard_singlet_S(params, 0.0, q)

"""
    hubbard_singlet_delta(params, K, q) -> Float64
    hubbard_singlet_delta(params, q)    -> Float64   (K = 0)

Exact singlet two-electron phase shift δ_K(q) = π/2 + atan(4 t cos(K/2) sin q / U).
"""
hubbard_singlet_delta(params::HubbardParams, K::Real, q::Real) =
    π/2 + atan(4 * params.t * cos(K/2) * sin(q) / params.U)
hubbard_singlet_delta(params::HubbardParams, q::Real) = hubbard_singlet_delta(params, 0.0, q)
