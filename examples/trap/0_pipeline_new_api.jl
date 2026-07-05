# =====================================================================
#  Trap pipeline on the library API (branch library-restructure).
#  Replaces the hand-rolled steps of 1–4: everything below is
#  model-independent library calls except the TFI block definition.
# =====================================================================
using Scattensor, SparseArrays, LinearAlgebra, Plots

# --- model block (the ONLY model-specific lines) ---------------------
d, L = 2, 12
J, h = 1.0, 1.5                                  # near-critical TFI
σx = sparse([0.0 1.0; 1.0 0.0]); σz = sparse([1.0 0.0; 0.0 -1.0])
id2 = sparse(1.0I, 2, 2)
H0 = -J * kron(σz, σz) - 0.5h * (kron(σx, id2) + kron(id2, σx))

# --- library pipeline ------------------------------------------------
H = summation_local(H0, d, L; pbc = true)
T = operator_translation(SparseMatrixCSC, d, L)

Hp = deformed_hamiltonian(H0, d, L; T = T)       # SSD trap  H' = Σ V(j) h_j
en, seeds = trap_seeds(Hp, 20)                   # 20 localized seeds, one Krylov run
H̃, T̃, _ = effective_pair(H, T, seeds, L)         # translated basis + Löwdin

eff = dispersion_relation(H̃, T̃, L; nlevels = 10) # effective bands (N·L problem)
exact = sector_spectrum(H, T, L; nlevels = 10)   # exact reference (fast sector ED)

E0eff = minimum(energy(s) for s in eff)
E0ex  = minimum(energy(s) for s in exact)
p = Plots.scatter([Float64(s.koverpi) for s in exact], [energy(s) - E0ex for s in exact];
    ms = 7, mc = :white, msc = :black, label = "exact (sector ED)",
    xlabel = "k/π", ylabel = "E − E0", dpi = 200,
    title = "Trap-effective vs exact, TFI h/J = $h, L = $L")
Plots.scatter!(p, [Float64(s.koverpi) for s in eff], [energy(s) - E0eff for s in eff];
    ms = 3, mc = :crimson, label = "trap effective (N = 20)")
mkpath(joinpath(@__DIR__, "plots"))
Plots.savefig(p, joinpath(@__DIR__, "plots", "pipeline_new_api.png"))
println("Done → plots/pipeline_new_api.png")
