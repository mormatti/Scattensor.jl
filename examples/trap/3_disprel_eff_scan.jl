
using SparseArrays
using Scattensor
using Revise
using Plots
using Colors
using LinearAlgebra
using KrylovKit
using Printf
using LaTeXStrings

# =====================================================================
#  Effective dispersion relation: cumulative basis scan N = 1..Nmax
# =====================================================================

# --- Parameters ---
d = 2
L = 12
L0 = 2
J = 1.0
h = 1.5

Nmax = 5
ic   = L ÷ 2
Vamp = 1.0

# --- Local operators ---
σx  = SparseMatrixCSC([0.0 1.0; 1.0 0.0])
σz  = SparseMatrixCSC([1.0 0.0; 0.0 -1.0])
id2 = SparseMatrixCSC([1.0 0.0; 0.0 1.0])
⊗(A, B) = kron(A, B)

H0 = -J * (σz ⊗ σz) - 0.5 * h * (σx ⊗ id2) - 0.5 * h * (id2 ⊗ σx)

# --- Build H, T once -------------------------------------------------
println("Constructing H, T ...")
H = summation_local(H0, d, L; pbc = true)
T = operator_translation(SparseMatrixCSC, d, L)

# --- Build H' with cosine trap ---------------------------------------
println("Constructing H' (cosine trap) ...")
V(i) = 1.0 - cos(2π * (i - ic) / L)
Idext = operator_identity(SparseMatrixCSC, d^(L - L0))
H0ext = kron(H0, Idext)
Hprime = spzeros(ComplexF64, d^L, d^L)
Hcur = SparseMatrixCSC{ComplexF64}(H0ext)
for j in 0:(L-1)
    global Hcur, Hprime
    Hprime += V(j) * Hcur
    Hcur = T * Hcur * T'
end
Hprime = (Hprime + Hprime') / 2

# --- Diagonalize H' for Nmax lowest states ---------------------------
println("Diagonalizing H' for $Nmax lowest levels ...")
en_p, st_p, _ = eigsolve(Hprime, size(Hprime, 1), Nmax, :SR, ComplexF64;
                         ishermitian = true,
                         krylovdim = max(2Nmax, Nmax + 20))
en_p = real.(en_p[1:Nmax])
st_p = st_p[1:Nmax]
println("Lowest $Nmax eigenvalues of H':")
for (i, e) in enumerate(en_p); @printf("  n=%2d   E_n' = %.8f\n", i, e); end

# =====================================================================
#  Helper: cumulative effective dispersion for given Ncur ≤ Nmax
# =====================================================================
function effective_drel(Ncur, st_p, H, T, L; tolfac = 1e-8)
    dimH = size(H, 1)
    Mdim = Ncur * L
    Vmat = zeros(ComplexF64, dimH, Mdim)
    for n in 1:Ncur
        Tpvn = Vector{ComplexF64}(st_p[n])
        for p in 0:(L-1)
            Vmat[:, (n - 1) * L + p + 1] = Tpvn
            Tpvn = T * Tpvn
        end
    end
    S   = Vmat' * Vmat;   S   = (S + S') / 2
    H_B = Vmat' * (H * Vmat); H_B = (H_B + H_B') / 2
    T_B = Vmat' * (T * Vmat)

    eS_vals, eS_vecs = eigen(Hermitian(S))
    eS_vals = real.(eS_vals)
    tol = tolfac * maximum(eS_vals)
    keep = eS_vals .> tol
    W = eS_vecs[:, keep] * Diagonal(1 ./ sqrt.(eS_vals[keep]))
    Htilde = W' * H_B * W;   Htilde = (Htilde + Htilde') / 2
    Ttilde = W' * T_B * W
    println("  N=$Ncur: Mdim=$Mdim, Meff=$(size(W,2)), ",
            "[H̃,T̃]=$(norm(Htilde*Ttilde - Ttilde*Htilde))")
    return dispersion_relation(Matrix(Htilde), Matrix(Ttilde), L; nlevels = 10)
end

# =====================================================================
#  Scan N = 1..Nmax, collect (k, E) per N
# =====================================================================
results = Vector{Tuple{Int, Vector{Float64}, Vector{Float64}}}()
for Ncur in 1:Nmax
    println()
    println("=== Effective basis with N = $Ncur seeds ===")
    drel = effective_drel(Ncur, st_p, H, T, L)
    ks  = Float64[]
    Es  = Float64[]
    for bs in drel
        # mirror to ±k for visualization, replicate to [-π, 2π]
        for kop in (bs.koverpi, -bs.koverpi)
            base = π * kop
            for shift in (-2π, 0.0, 2π)
                k = base + shift
                if -π - 1e-12 ≤ k ≤ 2π + 1e-12
                    push!(ks, k)
                    push!(Es, energy(bs))
                end
            end
        end
    end
    push!(results, (Ncur, ks, Es))
end

# =====================================================================
#  Plot: blue → red palette, larger markers for smaller N (drawn first)
# =====================================================================
palette = range(colorant"blue", stop = colorant"red", length = Nmax)
sizes   = range(9, stop = 4, length = Nmax)        # bigger for smaller N
alphas  = range(0.75, stop = 0.95, length = Nmax)

# Build axis limits from union of all data
all_E = vcat([r[3] for r in results]...)
Emin, Emax = extrema(all_E)
pad = 0.05 * (Emax - Emin)
ylims = (Emin - pad, Emax + pad)
xlims = (-π - 0.2, 2π + 0.2)

xticks_pos = [-π, -π/2, 0, π/2, π, 3π/2, 2π]
xticks_lab = [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]

pl = plot(;
    xlabel = L"k", ylabel = L"\mathcal{E}_k",
    title  = @sprintf("Effective spectrum, cumulative N (L=%d, h=%.2f)", L, h),
    legend = :topright, framestyle = :box,
    xlims = xlims, ylims = ylims,
    xticks = (xticks_pos, xticks_lab),
    size = (900, 650),
)

for (i, (Ncur, ks, Es)) in enumerate(results)
    scatter!(pl, ks, Es;
        color = palette[i], markersize = sizes[i], markerstrokewidth = 0.3,
        markerstrokecolor = :black, alpha = alphas[i],
        label = "N = $Ncur",
    )
end

plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)
output_path = joinpath(plots_dir, @sprintf("disprel_eff_scan_L%d_h%.2f.png", L, h))
savefig(pl, output_path)
println()
println("Plot saved to $output_path")
