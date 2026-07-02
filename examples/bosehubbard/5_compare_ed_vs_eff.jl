
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
#  BH (d=3) — Compare ED vs effective basis (N=2, 5, 10), for t=1 and t=2
# =====================================================================

# --- Shared parameters ---
d  = 3
L  = 10
L0 = 2
U  = 10.0
ε  = 100.0
ic = L ÷ 2
nlev = 20
Nseeds = 10

# --- Local operators ---
b   = SparseMatrixCSC(spdiagm(1 => [sqrt(i) for i in 1:(d-1)]))
bd  = SparseMatrixCSC(sparse(b'))
nop = SparseMatrixCSC(spdiagm(0 => Float64.(0:(d-1))))
id3 = SparseMatrixCSC(spdiagm(0 => ones(d)))
⊗(A, B) = kron(A, B)

# Build the BH 2-local H0 for given hopping t
function bh_H0(t)
    H0 = -t * (bd ⊗ b + b ⊗ bd)
    H0 += 0.5 * (U/2) * (nop * (nop - id3)) ⊗ id3
    H0 += 0.5 * (U/2) * id3 ⊗ (nop * (nop - id3))
    H0 += 0.5 * ε * (nop ⊗ id3 + id3 ⊗ nop)
    return H0
end

# Trap weight V(j) (cosine, min at ic)
Vtrap(j) = 1.0 - cos(2π * (j - ic) / L)

# Build the weighted H' = Σ_j Vtrap(j) · T^j H0ext T^{-j}
function build_Hprime(H0, T, L, L0)
    Idext = operator_identity(SparseMatrixCSC, d^(L - L0))
    H0ext = SparseMatrixCSC{ComplexF64}(kron(H0, Idext))
    Hprime = spzeros(ComplexF64, size(H0ext, 1), size(H0ext, 2))
    Hcur = H0ext
    for j in 0:(L-1)
        Hprime += Vtrap(j) * Hcur
        Hcur = T * Hcur * T'
    end
    return (Hprime + Hprime') / 2
end

# Effective dispersion from Ncur seeds (Löwdin with cutoff)
function effective_drel(Ncur, st_p, H, T, L; tolfac = 1e-8, nlev = 20)
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
    println("  Ncur=$Ncur  Mdim=$(Ncur*L)  Meff=$(size(W,2))  ",
            "[H̃,T̃]=$(norm(Htilde*Ttilde - Ttilde*Htilde))")
    return dispersion_relation(Matrix(Htilde), Matrix(Ttilde), L; nlevels = nlev)
end

# Expand a drel vector into plot points over [-π, 2π] (mirror to ±k)
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

# =====================================================================
#  Main loop over t values
# =====================================================================

plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)

xticks_pos = [-π, -π/2, 0, π/2, π, 3π/2, 2π]
xticks_lab = [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]

for t_val in (1.0, 2.0)
    println("\n##############################################")
    println("##  t = $t_val")
    println("##############################################")

    H0 = bh_H0(t_val)
    println("Building H, T ...")
    H = summation_local(H0, d, L; pbc = true)
    T = operator_translation(SparseMatrixCSC, d, L)

    println("Building H' (cosine trap) ...")
    Hprime = build_Hprime(H0, T, L, L0)

    println("Diagonalizing H' for $Nseeds lowest seeds ...")
    en_p, st_p, _ = eigsolve(Hprime, size(Hprime, 1), Nseeds, :SR, ComplexF64;
                             ishermitian = true,
                             krylovdim = max(2Nseeds, Nseeds + 25))
    st_p = st_p[1:Nseeds]
    println("Lowest $Nseeds eigenvalues of H':")
    for (i, e) in enumerate(real.(en_p[1:Nseeds]))
        @printf("  n=%2d   E_n' = %.4f\n", i, e)
    end

    println("\n=== Exact ED ===")
    drel_ed = dispersion_relation(H, T, L; nlevels = nlev)

    println("\n=== Effective N=2 ===")
    drel_N2 = effective_drel(2, st_p, H, T, L; nlev = nlev)
    println("\n=== Effective N=5 ===")
    drel_N5 = effective_drel(5, st_p, H, T, L; nlev = nlev)
    println("\n=== Effective N=10 ===")
    drel_N10 = effective_drel(10, st_p, H, T, L; nlev = nlev)

    k_ed,  E_ed  = expand_to_plot(drel_ed)
    k_N2,  E_N2  = expand_to_plot(drel_N2)
    k_N5,  E_N5  = expand_to_plot(drel_N5)
    k_N10, E_N10 = expand_to_plot(drel_N10)

    Emin, Emax = extrema(E_ed)
    pad = 0.05 * (Emax - Emin)
    ylims = (Emin - pad, Emax + pad)
    xlims = (-π - 0.2, 2π + 0.2)

    pl = plot(;
        xlabel = L"k", ylabel = L"\mathcal{E}_k",
        title = @sprintf("BH (d=3) ED vs eff.  L=%d, t=%.1f, U=%.1f, ε=%.1f", L, t_val, U, ε),
        legend = :topright, framestyle = :box,
        xlims = xlims, ylims = ylims,
        xticks = (xticks_pos, xticks_lab),
        size = (1000, 700),
    )
    scatter!(pl, k_ed, E_ed;
        color = :white, markerstrokecolor = :black, markerstrokewidth = 1.2,
        markersize = 11, label = "Exact ED")
    scatter!(pl, k_N2, E_N2;
        color = colorant"royalblue", markerstrokewidth = 0.0,
        markersize = 7, alpha = 0.85, label = "Eff. N=2")
    scatter!(pl, k_N5, E_N5;
        color = colorant"mediumpurple", markerstrokewidth = 0.0,
        markersize = 5, alpha = 0.90, label = "Eff. N=5")
    scatter!(pl, k_N10, E_N10;
        color = colorant"crimson", markerstrokewidth = 0.0,
        markersize = 3, alpha = 0.95, label = "Eff. N=10")

    fname = @sprintf("compare_bh_L%d_t%.1f_U%.1f_eps%.1f.png", L, t_val, U, ε)
    output_path = joinpath(plots_dir, fname)
    savefig(pl, output_path)
    println("Plot saved to $output_path")
end
