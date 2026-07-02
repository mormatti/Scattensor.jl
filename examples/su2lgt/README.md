# (1+1)D SU(2) lattice gauge theory — periodic ring, and why Lüscher stalls

This folder applies the Lüscher pipeline to the **spin/qubit formulation of the
(1+1)D SU(2) lattice gauge theory** of Barata–Hormaza–Kang–Qian
([arXiv:2511.00154](https://arxiv.org/abs/2511.00154), eq. 2.8), the model whose
real-time meson-meson scattering is shown in that paper's Figure 4.

## Model (periodic version)

`N` color-doubled sites → `2N` qubits (`φₙ = (ψ_{2n−1}, ψ_{2n})`, two colors).
After integrating out the gauge field, the Hamiltonian is

$$H = \frac{1}{2a}\sum_{j}\!\big(\sigma^+_j\sigma^-_{j+1}+\text{h.c.}\big)
   + ma\sum_n \tfrac{(-1)^n}{2}\big(\sigma^z_{2n-1}+\sigma^z_{2n}\big)
   + \frac{ag^2}{2}\sum_{n,m} G(n-m)\,\mathbf{Q}_n\!\cdot\!\mathbf{Q}_m ,$$

with color charges $Q_n^z=\tfrac14(\sigma^z_{2n-1}-\sigma^z_{2n})$,
$Q_n^+=\sigma^+_{2n-1}\sigma^-_{2n}$, and $G$ the **periodic** 1D-lattice Green's
function $(-\Delta)^{-1}$ — the translation-invariant analogue of the paper's
open-boundary linear (confining) Coulomb potential, so that lattice momentum is a
good quantum number. Baryon number $B=\tfrac14\sum_j\sigma^z_j$ is conserved; on a
ring only total-color-**singlet** states are physical, enforced with a
$\Lambda(\mathbf{Q}_\text{tot})^2$ penalty. Files: `su2_lgt.jl` (Hamiltonian),
`su2_scatter.jl` (momentum-resolved spectra by baryon sector). Fig-4 parameters:
`a=1, ga=5, ma=0.2` (`x=1/(ga)²=0.04`, deep strong coupling).

## What we validated

* Vacuum energy `E0 = −maN` (staggered Dirac sea), correct color charges.
* Lowest **B=0** excitation = **baryon–antibaryon** (doublon–holon), `ΔE = 4m`,
  color-neutral (`Σ⟨Q·Q⟩ = 0`) — *not* the meson.
* The **meson** (q-q̄, `Σ⟨Q·Q⟩ = 3/2`) sits at `ΔE ≈ 2m + 3ag²/8 ≈ 8`, a heavy
  resonance buried in the lighter baryon-antibaryon continuum.
* The single **baryon** (B=1, the lightest *particle*, `ΔE = 2m`) has an
  **exactly flat band**: bandwidth `≈ 4×10⁻¹³` at every coupling.

## Why the Lüscher phase shift can't be extracted here

Lüscher inverts the finite-volume relation between a **two-particle energy** and
the **relative momentum** `q`, via `qL + 2δ(q) = 2πm` with
`E_free(K,q) = ε(K/2+q) + ε(K/2−q)`. This needs the one-particle dispersion
`ε(k)` to **vary with `k`**. In this theory:

* the lightest particle (baryon) is an **exactly localized** compact state —
  `ε(k) = const`, so `E_free` carries no information about `q`; the quantization
  condition degenerates and no phase shift exists to extract;
* the meson *does* carry energy but is heavy and embedded in a lighter continuum
  (it can decay), so it is not a clean asymptotic state for energy-level Lüscher.

Physically: the hadrons of this confining lattice theory are **heavy / static**.
That is precisely why the paper scatters them with **real-time wavepackets**
(TEBD) — a method that imparts momentum to even a static particle by construction
— rather than from the low-energy spectrum. The energy-based Lüscher method that
works cleanly for the dispersing magnons/electrons of the spin chains
(`../luscher/`: XXZ, Hubbard, Ising) does **not** transfer to these hadrons.

## Natural next step

To reproduce Figure 4's physics one would use the **real-time** route already
half-present in this repo: prepare baryon/meson wavepackets, time-evolve, and read
the **Wigner time delay** `τ = 2 dδ/dE` from the collision — the dynamical
counterpart of the phase shift, valid for static/heavy particles.
