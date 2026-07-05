# Library restructure — inventory and promotion map

Branch `library-restructure`, July 2026. Rule for promotion into `src/`:
**model-free signature** (takes `(H0, d, L)`, `(H, T)`, `(ψ, Ω)` — never a
model name) **+ used ≥ 2×** across projects **+ testable on a solvable model**.
Everything hardcoded (couplings, file paths, plot cosmetics) stays in
`examples/` or in the downstream projects (`digitization`, `spectroscopy`).

## Promoted in this branch

| module (src/scattensor_functions/) | contents | origin |
|---|---|---|
| `sectors.jl` | monomial-symmetry block decomposition (orbits → Bloch isometries), `momentum_basis`, `sector_hamiltonian`, `sector_spectrum` (fast exact replacement for full-space projection in `dispersion_relation`) | new (request 2026-07-05) |
| `localizability.jl` | reduced transition operator trace norm (QR trick), `localizability` C̃(ℓ), `support_size` ℓ_loc, `optimal_creator` φ* | digitization scripts 36/38/39 |
| `trap.jl` | `ssd_weight`, `deformed_hamiltonian` (Σ V(j)h_j), `trap_seeds`, `effective_pair` (translated basis + Löwdin) | examples/trap 1–4 |
| `fits.jl` | `wrap_mod_pi`, profile-linear `fit_arctan` | digitization scripts 27/28/32/40 |
| `tdvp_time_evolution.jl` (upgraded) | `save_states` switch, `plots_every`, atomic checkpoints + `tdvp_checkpoint`, live `dashboard.html` | pain points of July 2026 runs |

## Known duplications in downstream scripts (to be replaced by library calls)

- `embed(ops, j0)` small-operator kron placement: digitization 22, 34, 36 (×3+)
  → candidate `insert_local`-style matrix helper (TODO: unify with the
  existing ITensor `insert_local`).
- `packet_op(jc, k)` Gaussian-packet creation via `summation_local` with
  convolution: digitization 26, 30, 33, 37 → thin helper over the existing
  library function (TODO).
- per-script arctan grids (27/28/40) → `fit_arctan`.
- `tracenorm`/`locC`/`ell_loc` (36/38/39) → `localizability.jl`.
- cosine-trap construction (examples/trap 1–4) → `deformed_hamiltonian`.

## Not promoted (and why)

- DP/Viterbi phase unwrap (script 31/32): flattens genuine π rises — kept as
  a script-level experiment, not a library contract.
- Prony toolkit (script 23): promising but single-use so far; revisit.
- Model definitions (`su2_ladder_block`, XXZ, …): live in the projects'
  `src/models.jl`, never in the library.

## TODO next on this branch

- [ ] `examples/trap` rewritten on top of `deformed_hamiltonian`/`effective_pair`.
- [ ] port digitization 34–40 to `spectroscopy/scripts` as thin drivers.
- [ ] product-symmetry sectors (translation × Z₂) — compose `symmetry_blocks`.
- [ ] `sector_spectrum` for the `dispersion_relation` TN (MPS) path.
- [ ] checkpoint auto-resume helper (`tdvp_time_evolution!(…, resume = path)`).
