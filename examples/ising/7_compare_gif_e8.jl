
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
#  TFI + longitudinal field near critical point (E8 region)
#  H = -J Σ σz σz - h Σ σx - hz Σ σz
#  At J=h=1 the TFI is critical; a small hz drives the system to the
#  E8 (Zamolodchikov) integrable point with 8 distinct massive particles.
#
#  Same machinery as BH gif: trap + Wannier-translated effective basis,
#  Löwdin orthonormalization, dispersion_relation, variance coloring.
# =====================================================================

d  = 2
L  = 13
L0 = 2
J  = 1.0
h  = 1.0    # TFI at critical point
hz = 0.3    # longitudinal field, more distant from critical point
ic = L ÷ 2
nlev = 40
Nseeds = 50

σx  = SparseMatrixCSC([0.0 1.0; 1.0 0.0])
σz  = SparseMatrixCSC([1.0 0.0; 0.0 -1.0])
id2 = SparseMatrixCSC([1.0 0.0; 0.0 1.0])
⊗(A, B) = kron(A, B)

H0  = -J * (σz ⊗ σz)
H0 += -0.5 * h  * (σx ⊗ id2 + id2 ⊗ σx)
H0 += -0.5 * hz * (σz ⊗ id2 + id2 ⊗ σz)

Vtrap(j) = 1.0 - cos(2π * (j - ic) / L)

function build_Hprime(H0, T, L, L0)
    Idext = operator_identity(SparseMatrixCSC, d^(L - L0))
    H0ext = SparseMatrixCSC{ComplexF64}(kron(H0, Idext))
    Hp = spzeros(ComplexF64, size(H0ext, 1), size(H0ext, 2))
    Hcur = H0ext
    for j in 0:(L-1)
        Hp += Vtrap(j) * Hcur
        Hcur = T * Hcur * T'
    end
    return (Hp + Hp') / 2
end

function effective_drel(Ncur, st_p, H, T, L; tolfac = 1e-8, nlev = 40)
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
    drel = dispersion_relation(Matrix(Htilde), Matrix(Ttilde), L; nlevels = nlev)

    VW = Vmat * W
    variances = Float64[]
    for bs in drel
        ψ = VW * bs.data
        Hψ = H * ψ
        eH  = real(dot(ψ, Hψ))
        eH2 = real(dot(Hψ, Hψ))
        push!(variances, eH2 - eH^2)
    end
    return drel, variances
end

function expand_to_plot(drel, vals = nothing)
    ks, Es = Float64[], Float64[]
    cs = vals === nothing ? nothing : Float64[]
    for (i, bs) in enumerate(drel)
        for kop in (bs.koverpi, -bs.koverpi)
            base = π * kop
            for shift in (-2π, 0.0, 2π)
                k = base + shift
                if -π - 1e-12 ≤ k ≤ 2π + 1e-12
                    push!(ks, k); push!(Es, energy(bs))
                    vals === nothing || push!(cs, vals[i])
                end
            end
        end
    end
    return vals === nothing ? (ks, Es) : (ks, Es, cs)
end

plots_dir = joinpath(@__DIR__, "plots")
mkpath(plots_dir)

xticks_pos = [-π, -π/2, 0, π/2, π, 3π/2, 2π]
xticks_lab = [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]

println("Building H, T ...")
H = summation_local(H0, d, L; pbc = true)
T = operator_translation(SparseMatrixCSC, d, L)

println("Building H' (cosine trap) ...")
Hp = build_Hprime(H0, T, L, L0)

println("Diagonalizing H' for $Nseeds lowest seeds ...")
en_p, st_p, _ = eigsolve(Hp, size(Hp, 1), Nseeds, :SR, ComplexF64;
                         ishermitian = true,
                         krylovdim = max(2Nseeds, Nseeds + 30))
st_p = st_p[1:Nseeds]

println("Computing exact ED dispersion ...")
drel_ed = dispersion_relation(H, T, L; nlevels = nlev)
k_ed, E_ed = expand_to_plot(drel_ed)

Emin, Emax = extrema(E_ed)
pad = 0.05 * (Emax - Emin)
ylims = (Emin - pad, Emax + pad)
xlims = (-π - 0.2, 2π + 0.2)

clims = (-12.0, 2.0)

println("Building animation frames ...")
anim = Animation()
for Ncur in 1:Nseeds
    print("  frame N=$Ncur ... ")
    drel_N, vars_N = effective_drel(Ncur, st_p, H, T, L; nlev = nlev)
    logvars = log10.(abs.(vars_N) .+ 1e-30)
    k_N, E_N, c_N = expand_to_plot(drel_N, logvars)
    println("done")

    pl = plot(;
        xlabel = L"k", ylabel = L"\mathcal{E}_k",
        title = @sprintf("N = %d", Ncur),
        legend = :topright, framestyle = :box,
        xlims = xlims, ylims = ylims,
        xticks = (xticks_pos, xticks_lab),
        size = (500, 1000),
        colorbar_title = L"\log_{10}|\mathrm{var}(H)|",
    )
    scatter!(pl, k_ed, E_ed;
        color = :black, markerstrokewidth = 0.0,
        markersize = 3, label = "Exact ED")
    scatter!(pl, k_N, E_N;
        zcolor = c_N, clims = clims, c = :plasma,
        markershape = :xcross,
        markersize = 7, markerstrokewidth = 2.0,
        label = "Eff. N=$Ncur")

    frame(anim, pl)
end

fname = @sprintf("compare_ising_e8_gif_L%d_J%.1f_h%.1f_hz%.2f.gif", L, J, h, hz)
output_path = joinpath(plots_dir, fname)
gif(anim, output_path, fps = 3)
println("GIF saved to $output_path")
