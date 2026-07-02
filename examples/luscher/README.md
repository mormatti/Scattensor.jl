# Lüscher phase-shift extraction for the periodic Ising chain

A modular Julia pipeline that extracts a nonperturbative **two-particle elastic
scattering phase shift** `δ_K(q)` from finite periodic spin chains using a
**1D Lüscher-style method**, built on top of the `Scattensor` library.

The pipeline is model-agnostic. Three models ship with it:

```
Ising:    H = -J Σ Z_j Z_{j+1} - h_x Σ X_j - h_z Σ Z_j            (d=2, the assignment model)
XXZ:      H = -J Σ (SˣSˣ + SʸSʸ + Δ SᶻSᶻ) - h Σ Sᶻ                (d=2, analytic benchmark)
Hubbard:  H = -t Σ c†c + U Σ n↑n↓ + μ Σ n + (isolation penalties)  (d=4, analytic benchmark)
```

with PBC. The XXZ and Hubbard chains are Bethe-ansatz integrable and provide the
cleanest tests: their two-particle S-matrices are known in closed form (see
Benchmarks). Adding a model only requires a `local_hamiltonian_term(params)`
method (and, for non-product channels, a parity operator).

## Physics chain

| step | quantity | formula |
|------|----------|---------|
| 1 | finite-volume spectrum | diagonalize `H`, label by `K = 2π Kint/L` |
| 2 | vacuum | `E0(L)` = lowest level |
| 3 | one-particle band | `ε_L(k) = E_1(L,k) − E_0(L)` |
| 4 | dispersion fit | `ε(k) = c_0 + Σ_r c_r cos(r k)` |
| 5 | free two-particle | `E_free(K,q) = ε(K/2+q) + ε(K/2−q)`, `q_m = 2π m/L` |
| 6 | relative momentum | solve `ΔE = ε(K/2+q) + ε(K/2−q)` |
| 7 | phase shift | `δ_K(q) = π m − qL/2`, `S_K(q) = e^{2iδ}` |
| 8 | global fit | `δ_K(q;θ) = a_0 + a_1 q + a_2 q² + a_3 q³`, minimize `D(θ)` |

No operator overlaps are used: bands and two-particle levels are identified from
energy–momentum structure, continuity in `K`, level ordering, and the
two-particle threshold only.

## Files

| file | contents |
|------|----------|
| `Model.jl` | `IsingParams`, `XXZParams`, `LuscherConfig`, **library-isolation wrappers** (`build_hamiltonian`, `build_translation`, `raw_dispersion`) |
| `Spectrum.jl` | `Level`, `compute_momentum_spectrum`, vacuum identification |
| `Dispersion.jl` | one-particle band (with false-vacuum guard), `fit_dispersion_fourier`, `epsilon` |
| `Luscher.jl` | matching, `solve_relative_momentum`, `extract_phase_shift`, unwrapping, global fit, scans |
| `Benchmark.jl` | exact references: free-fermion dispersion, E8 masses + `S₁₁(θ)`, **XXZ two-magnon S-matrix** |
| `Plotting.jl` | all figures |
| `run_luscher_ising.jl` | end-to-end driver (the four scenarios below) |

The wrappers in `Model.jl` are the **only** place the `Scattensor` library is
touched — re-point them to swap the diagonalization backend.  The pipeline is
model-agnostic; a new model only needs a `local_hamiltonian_term(params)` method.

## Run

```bash
julia --project=. examples/luscher/run_luscher_ising.jl
```

Outputs are organized per scenario: `plots/<scenario>/*.png` and
`results/<scenario>/{phase_points.csv, summary.json}`.

## Benchmarks (the point of the exercise)

### ★ XXZ easy-axis ferromagnet — fully analytic δ(q)  *(headline)*

The XXZ chain `H = -J Σ (SˣSˣ + SʸSʸ + Δ SᶻSᶻ) - h Σ Sᶻ` is **Bethe-ansatz
integrable for every Δ**.  In the easy-axis ferromagnet (Δ>1) the vacuum is the
exact product state |↑↑…↑⟩, single magnons have ε(k)=J(Δ−cos k)+h, and the
two-magnon elastic S-matrix at total momentum K=0 is known in **closed form**:

```
S(q) = − (1 − Δ e^{−iq}) / (1 − Δ e^{iq}),   δ(q) = π/2 + atan2(Δ sin q, 1 − Δ cos q).
```

The uniform field `h` shifts an n-magnon state by `n·h` (it lifts the
two-magnon bound states / continuum above the one-magnon band so the dispersion
is cleanly identified) **without changing the scattering phase**.  The extracted
`S(q)` matches the exact formula to ~machine precision on the clean levels —
this is the real validation of the Lüscher extraction.  (At Δ=0 it reduces to
free fermions, `S=−1`.)  Ref: Bethe (1931); Karbach–Müller, cond-mat/9809162.

### ★ 1D Hubbard — fully analytic two-electron δ(q)

The Hubbard chain (Lieb–Wu integrable) is benchmarked in the **two-electron
sector over the empty vacuum**: one electron is free, ε(k)=−2t cos k+μ, and an
opposite-spin pair scatters via U.  The interacting (spin-singlet) two-body
S-matrix at K=0 is **closed form**:

```
δ(q) = π/2 + atan(4 t sin q / U),   S(q) = −(U + 4 i t sin q)/(U − 4 i t sin q).
```

Limits: U→0 ⇒ S=+1 (free); U→∞ ⇒ S=−1 (fermionized hard core). Ref: Lieb & Wu,
PRL 20, 1445 (1968); Essler et al., *The 1D Hubbard Model*.

Isolating this channel needs care (and is a nice illustration of the method):
a chemical potential μ makes the empty state the vacuum; cheap diagonal penalties
`(Sᶻ_tot)² P_{N=2}` and `P_{N≥3}` remove the 2↑/2↓ and ≥3-particle sectors; and
the **spin singlet (sees U) vs triplet (free)** are separated by **reflection
parity at K=0** (the hard-core-boson representation used here carries no SU(2),
so reflection — not S²_tot — is the right separator).  With that, the extracted
`S(q)` matches the exact formula to ~1e-11.  d=4 limits ED to small `L` (6–9).

### Wigner time delay

For every model with an exact δ, the pipeline also computes the **Wigner time
delay** τ(q) = 2 dδ/dE (E the two-particle energy): the exact curve from the
analytic δ(q), and the data estimate from a discrete central derivative of the
extracted points (`timedelay.png`). Free particles give τ=0; repulsive
interactions give a threshold-peaked τ(q).

### Nonzero total momentum (multi-K)

By default each scenario uses the centre-of-mass sector `K = 0`. The pipeline is
**K-general**, and XXZ, the free fermion and the disordered point also extract
nonzero total-momentum sectors `K = 2π Kint/L` (`phaseshift_multiK.png`). Two
subtleties make `K ≠ 0` non-trivial:

* **Bethe parity shift.** Two free particles have `k₁,₂ = 2π j₁,₂/L` with
  `j₁+j₂ = Kint`, so the relative momentum `q = π(j₁−j₂)/L` carries the *parity of
  `Kint`*. For **even** `Kint` the free relative momenta are `2πm/L`, but for **odd**
  `Kint` they shift to `π(2m−1)/L` (the moving-frame / Neveu–Schwarz half-integer
  shift). The quantization becomes `δ_K(q) = πm − π·(Kint mod 2)/2 − qL/2`.
* **K-dependent S-matrix.** The two-body phase depends on `K`. For XXZ it factorizes
  cleanly, `S_K(q) = −(cos(K/2) − Δe^{−iq})/(cos(K/2) − Δe^{iq})`, and for Hubbard the
  coupling enters as `sin k₁ − sin k₂ = 2cos(K/2) sin q`. Because `K = 2π Kint/L`
  depends on `L`, the multi-K plot is made at a single `L` (the largest not divisible
  by 4, to avoid the `q=π/2` artifact). The extracted points land on the K-dependent
  exact `δ_K(q)` to ~machine precision on the clean low levels — including the **odd**
  sectors, validating the parity shift.

Hubbard stays `K = 0`: its singlet/triplet separation uses **reflection parity**,
which is only a good quantum number at `K = 0` (reflection maps `K → −K`).

### Ising benchmarks

1. **Free-fermion point** `h_x=4, h_z=0` (deep in the disordered phase). Exact dispersion
   `ε(k)=2√(J²+h_x²−2Jh_x cos k)`, constant `S=−1` ⇒ `δ=π/2`. The gap `2(h_x−1)=6`
   exceeds the bandwidth (≈4), so the one-particle band is well separated from the
   4-particle continuum and the low levels are clean. (Near the critical point,
   `h_x≈1.5`, the bands overlap: the TFIM is a *paired* BdG free theory with no
   particle-number conservation, so 4-/6-quasiparticle states sit just above the
   2-particle threshold and contaminate the single-channel extraction — a textbook
   illustration of the method's inelastic-contamination warning. Also, the √-dispersion
   is not a finite cosine series, unlike XXZ/Hubbard, so its Fourier fit is only
   approximate.) Validates bookkeeping.
2. **Disordered + field** `h_x=4, h_z=1` (paramagnet, longitudinal field on). This is
   the free-fermion point with `h_z` switched on: still deep in the disordered phase,
   so the one-particle band stays **fully separated** (gap ≈ 6.8 > bandwidth ≈ 3.6, band
   top 10.4 < two-particle threshold 13.7), but the magnons now genuinely **interact**.
   The model is **non-integrable** (no closed-form S-matrix), yet because the band is
   clean the Lüscher extraction is well-behaved: the extracted `δ_K(q)` is a smooth curve
   that deviates from the free `π/2` by the interaction shift (≈0.2–0.3 rad at `h_z=1`).
   This is the contrast with the E8 point — *away* from criticality the band separates
   and the method works cleanly even without integrability. Kept for the clean-extraction
   demonstration; `δ` deviation from `π/2` is the physical content.
3. **E8 point** `h_x=1, h_z=0.04` (critical TFIM + field): Zamolodchikov's E8 theory.
   The **mass ratio m₂/m₁ → 2cos(π/5)=1.6180…** (golden ratio) converges cleanly
   with `L` (≈0.1 % at `L=14`). Refs: Zamolodchikov 1989; Coldea et al., *Science*
   **327**, 177 (2010). The exact `S₁₁(θ)={2/3}{2/5}{1/15}` is also overlaid, but a
   quantitative phase match needs larger `L`/smaller `h_z` than ED at `L≲14` allows.
   - **E8 isolated** `h_x=1, h_z=0.5` (scenario `e8_isolated`). The E8 *field theory*
     is integrable at any field, but on the lattice E8 is the small-`h_z` *scaling
     limit*. Raising `h_z` makes the gap `~O(1)` lattice scale, which **isolates the
     one-particle band** (band top ≈5.4 < `2m₁`≈7.4) and yields a *clean* phase
     extraction that tracks the E8 `δ₁₁(θ)` curve — but at the price of leaving the
     scaling limit: the golden ratio drifts to ≈1.4 % and the phase to ~0.1–1 rad.
     A crossover ("approximately E8 with an isolated band"), complementary to the
     pristine-but-contaminated `e8` point. The band-isolation ⇄ E8-precision tension
     is intrinsic: isolation needs gap ~ cutoff, E8 needs gap ≪ cutoff.
4. **Confinement point** `h_x=0.6, h_z=0.05` (assignment defaults). Ordered phase:
   kinks confine into mesons (McCoy–Wu), **non-integrable**, no exact `δ(q)`. The
   second ferromagnetic vacuum sits ≈`2h_zL⟨Z⟩` above the ground state; the band
   builder detects it and selects the genuine meson-at-rest by continuity, with
   warnings. Kept for physics, not validation.

### Why not XYZ?

The XYZ chain `H = Σ (Jₓ SˣSˣ + J_y SʸSʸ + J_z SᶻSᶻ)` is integrable (Baxter /
eight-vertex), but it is a **poor ED-Lüscher benchmark**: (i) its exact S-matrix
is **elliptic** (Jacobi theta functions); (ii) for Jₓ≠J_y the SˣSˣ−SʸSʸ terms do
not conserve Sᶻ, so there is **no product vacuum** and no simple magnon sector;
(iii) the elementary excitations are **topological kinks** (like Ising
confinement), so the one-/two-particle identification breaks at small `L`. Its
two clean limits are exactly the ones already covered here — Jₓ=J_y is XXZ, and
J_z=0 is the free-fermion XY model (`S=−1`). It is therefore omitted by design.
