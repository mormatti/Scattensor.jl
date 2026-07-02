
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
#  Compare three spectra: exact ED, effective N=2, effective N=5
# =====================================================================

# --- Parameters ---
d = 2
L = 12
L0 = 2
J = 1.0
h = 1.5
ic = L ÷ 2
Vamp = 1.0
nlev = 10

# --- Local operators ---
σx  = SparseMatrixCSC([0.0 1.0; 1.0 0.0])
σz  = SparseMatrixCSC([1.0 0.0; 0.0 -1.0])
id2 = SparseMatrixCSC([1.0 0.0; 0.0 1.0])
⊗(A, B) = kron(A, B)

H0 = -J * (σz ⊗ σz) - 0.5 * h * (σx ⊗ id2) - 0.5 * h * (id2 ⊗ σx)

println("Constructing H, T ...")
H = summation_local(H0, d, L; pbc = true)
T = operator_translation(SparseMatrixCSC, d, L)

# --- H' with cosine trap ---
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

# --- Diagonalize H' for 5 lowest seeds ---
Nseeds = 20
println("Diagonalizing H' for $Nseeds lowest levels ...")
en_p, st_p, _ = eigsolve(Hprime, size(Hprime, 1), Nseeds, :SR, ComplexF64;
                         ishermitian = true,
                         krylovdim = max(2Nseeds, Nseeds + 20))
st_p = st_p[1:Nseeds]

# --- Effective dispersion helper ---
function effective_drel(Ncur, st_p, H, T, L; tolfac = 1e-8, nlev = 10)
    dimH = size(H, 1)
    Vmat = zeros(ComplexF64, dimH, Ncur * L)
    for n in 1:Ncur
        Tpvn = Vector{ComplexF64}(st_p[n])
        for p in 0:(L-1)
            Vmat[:, (n - 1) * L + p + 1] = Tpvn
            Tpvn = T * Tpvn
        end
    end
    S   = Vmat' * Vmat;        S   = (S + S') / 2
    H_B = Vmat' * (H * Vmat);  H_B = (H_B + H_B') / 2
    T_B = Vmat' * (T * Vmat)
    eS_vals, eS_vecs = eigen(Hermitian(S))
    eS_vals = real.(eS_vals)
    keep = eS_vals .> tolfac * maximum(eS_vals)
    W = eS_vecs[:, keep] * Diagonal(1 ./ sqrt.(eS_vals[keep]))
    Htilde = W' * H_B * W; Htilde = (Htilde + Htilde') / 2
    Ttilde = W' * T_B * W
    return dispersion_relation(Matrix(Htilde), Matrix(Ttilde), L; nlevels = nlev)
end

# --- Exact ED dispersion ---
println()
println("=== Exact ED ===")
drel_ed = dispersion_relation(H, T, L; nlevels = nlev)

println()
println("=== Effective N=10 ===")
drel_N2 = effective_drel(10, st_p, H, T, L; nlev = nlev)

println()
println("=== Effective N=20 ===")
drel_N5 = effective_drel(20, st_p, H, T, L; nlev = nlev)

# --- Build (k, E) lists with ±k mirror and replicas in [-π, 2π] ---
function expand_to_plot(drel)
    ks, Es = Float64[], Float64[]
    for bs in drel
        for kop in (bs.koverpi, -bs.koverpi)
            base = π * kop
            for shift in (-2π, 0.0, 2π)
                k = base + shift
                if -π - 1e-12 ≤ k ≤ 2π + 1e-12
                    push!(ks, k); push!(Es, energy(bs))
                end
            end
        end
    end
    return ks, Es
end

k_ed, E_ed = expand_to_plot(drel_ed)
k_N2, E_N2 = expand_to_plot(drel_N2)
k_N5, E_N5 = expand_to_plot(drel_N5)

# Y-range from the exact set (so high spurious eff levels don't blow the axis)
Emin, Emax = extrema(E_ed)
pad = 0.05 * (Emax - Emin)
ylims = (Emin - pad, Emax + pad)
xlims = (-π - 0.2, 2π + 0.2)

xticks_pos = [-π, -π/2, 0, π/2, π, 3π/2, 2π]
xticks_lab = [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]

pl = plot(;
    xlabel = L"k", ylabel = L"\mathcal{E}_k",
    title  = @sprintf("ED vs effective (L=%d, h=%.2f)", L, h),
    legend = :topright, framestyle = :box,
    xlims = xlims, ylims = ylims,
    xticks = (xticks_pos, xticks_lab),
    size = (1000, 700),
)

# Exact ED: large hollow black circles (background reference)
scatter!(pl, k_ed, E_ed;
    color = :white, markerstrokecolor = :black, markerstrokewidth = 1.4,
    markersize = 11, label = "Exact ED")
# Effective N=2: medium blue filled
scatter!(pl, k_N2, E_N2;
    color = colorant"royalblue", markerstrokewidth = 0.0,
    markersize = 6, alpha = 0.85, label = "Eff. N=10")
# Effective N=20: small red filled (drawn on top)
scatter!(pl, k_N5, E_N5;
    color = colorant"crimson", markerstrokewidth = 0.0,
    markersize = 3, alpha = 0.95, label = "Eff. N=20")

plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)
output_path = joinpath(plots_dir,
    @sprintf("disprel_compare_L%d_h%.2f.png", L, h))
savefig(pl, output_path)
println()
println("Plot saved to $output_path")
