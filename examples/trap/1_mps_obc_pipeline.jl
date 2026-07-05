# =====================================================================
#  MPS/OBC production pipeline of the bulk-reconstruction spectroscopy
#  (branch tn-port). The full method — trap seeds, transition-operator
#  creators, bulk couplings, Wannier bands — with NO exact diagonalization
#  anywhere: plain finite DMRG + local operators + linear algebra, on an
#  OPEN chain. This is the mode that scales (L = 24 here; nothing below
#  is exponential in L).
#  Benchmark: transverse-field Ising, ε(k) = 2√(J² + h² − 2Jh cos k).
# =====================================================================
using Scattensor, SparseArrays, LinearAlgebra, Random, Plots
using ITensors, ITensorMPS
gr()
Random.seed!(7)

const d = 2
const J, h = 1.0, 1.5
const L = 24
σx = sparse([0.0 1.0; 1.0 0.0]); σz = sparse([1.0 0.0; 0.0 -1.0])
id2 = sparse(1.0I, 2, 2)
H0 = -J * kron(σz, σz) - 0.5h * (kron(σx, id2) + kron(id2, σx))
eps_exact(k) = 2 * sqrt(J^2 + h^2 - 2J * h * cos(k))

# ---------- DMRG vacuum on the open chain ----------
t0 = time()
Hmpo = summation_local(mpo_from_matrix(Matrix(H0), d), L)
sites = siteinds_main(Hmpo)
E0, Ω = dmrg(Hmpo, random_mps(sites; linkdims = 8);
             nsweeps = 12, maxdim = [10, 20, 40, 80, 120], cutoff = 1e-11,
             outputlevel = 0)
println("DMRG vacuum: E0 = $(round(E0, digits = 8))  ($(round(time() - t0, digits = 1)) s)")

# ---------- parabolic trap + DMRG seeds ----------
t1 = time()
trap = deformed_hamiltonian_mpo(H0, d, L; weight = parabolic_weight(L))
replace_siteinds!(trap, sites)
en, seeds = dmrg_seeds(trap, sites, 8; nsweeps = 12, outputlevel = 0)
ord = sortperm(en); en = en[ord]; seeds = seeds[ord]   # penalty-DMRG can return out of order
println("trap seeds: ", round.(en, digits = 4), "  ($(round(time() - t1, digits = 1)) s)")

# ---------- transition-operator creators (window scan) ----------
const ℓφ = 5
function seed_creator(seed)
    a = inner(Ω, seed)
    sP = add(seed, -a * Ω; cutoff = 1e-12)
    normalize!(sP)
    span = 4:(L - ℓφ - 2)
    tns = [sum(svdvals(mps_transition_operator(sP, Ω, s:(s + ℓφ - 1)))) for s in span]
    sbest = span[argmax(tns)]
    println("  creator: window at site $sbest, ‖A‖₁ = $(round(maximum(tns), digits = 4))")
    A = mps_transition_operator(sP, Ω, sbest:(sbest + ℓφ - 1))
    A ./ norm(A)
end
t2 = time()
creators = [seed_creator(seeds[n]) for n in 2:8]
println("creators extracted  ($(round(time() - t2, digits = 1)) s)")

# ---------- bulk couplings + Wannier interpolation ----------
t3 = time()
S, Hc = mps_bulk_couplings(Hmpo, Ω, creators, L ÷ 2, 6; cutoff = 1e-12, maxdim = 160)
println("bulk couplings  ($(round(time() - t3, digits = 1)) s)")
println("|S(r)| decay (α=1): ", round.(abs.(S[1, 1, 7:end]), sigdigits = 2))
ks, bands = wannier_bands(S, Hc; nk = 400, gram_tol = 1e-3)

# ---------- compare with the exact band ----------
kk = range(0, π; length = 300)
devs = Float64[]
for k in kk
    i = argmin(abs.(ks .- k))
    isnan(bands[i, 1]) && continue
    push!(devs, abs(bands[i, 1] - eps_exact(k)) / eps_exact(k))
end
println("max relative deviation from exact band: ", round(maximum(devs), sigdigits = 3))

p = Plots.plot(kk ./ π, eps_exact.(kk); lw = 2.5, color = :black, ls = :dash,
    label = "exact free fermions", xlabel = "k/π", ylabel = "E − E0",
    title = "ITF band from MPS/OBC bulk reconstruction, L = $L", dpi = 200,
    legend = :topleft)
sel = ks .>= 0
Plots.plot!(p, ks[sel] ./ π, bands[sel, 1]; lw = 2, color = :crimson,
    label = "MPS operator frame (OBC, no ED)")
mkpath(joinpath(@__DIR__, "plots"))
Plots.savefig(p, joinpath(@__DIR__, "plots", "mps_obc_pipeline.png"))
println("Done → plots/mps_obc_pipeline.png   (total $(round(time() - t0, digits = 1)) s)")
