# Tensor Networks

""" Computes the entanglement entropy in a position j for a ITensor MPS.
    """
function entanglement_entropy(ψ::MPS, j::Int)::Float64
    orthogonalize!(ψ, j)
    _,S = svd(ψ[j], (linkind(ψ,j), siteind(ψ,j+1)))
    return -sum(p^2 * log(p^2) for p in diag(S)) / log(2)
end

""" Computes the entanglement entropy array for a ITensor MPS.
    """
function entanglement_entropy(ψ::MPS)::Vector{Float64}
    L = length(ψ)
    Sᵥ(j) = entanglement_entropy(ψ, j)
    return [Sᵥ(j) for j in 1:L-1]
end

""" Apply the reflection of a finite MPS with uniform local dimension.
    """
function reflect(ψ::MPS)
    N, Sd = length(ψ), siteinds(ψ)
    ϕ = MPS(N)
    for j in 1:N
        h = N - j + 1
        ϕ[h] = ψ[j] * delta(Sd[j], Sd[h]')
    end
    noprime!(siteinds, ϕ)
    return ϕ
end

""" Apply the translation of a finite MPS with uniform local dimension.
    The translation is performed swapping couple of physical indices consecutively.
    """
function translate(ψ::MPS; dir = "right", cutoff = 1e-10)
    L = length(ψ)
    ϕ = copy(ψ)
    if dir == "left"
        for j in 1:L-1
            ϕ = swapbondsites(ϕ, j, cutoff = cutoff)
        end
    elseif dir == "right"
        for j in L-1:-1:1
            ϕ = swapbondsites(ϕ, j, cutoff = cutoff)
        end
    else
        error("Invalid direction")
    end
    return ϕ
end

# To polish
""" Computes the time evolution of a MPS with a MPO Hamiltonian.
    """
function tdvp_time_evolution(
    H::MPO,
    ψ::MPS,
    ψ₀::MPS,
    dt::Float64,
    Δt::Float64
    )

    N = Integer(Δt / dt) # Number of steps

    # We create the list of times 
    times = [n * dt for n in 1:N]
    positions = [x for x in 2:L-1]

    # We create a 2D array E of dimensions (N, L-2)
    En = zeros(N, L-2)
    EnLog = zeros(N, L-2)
    Linkdims = zeros(N, L-2)
    Maxlinkdim = zeros(N)

    for n in 1:N
        ψ = tdvp(H, ψ, -im * dt, maxlinkdim = 100, cutoff = 1e-10)
        normalize!(ψ)

        timeElapsed = @elapsed begin
        for j in 2:(L-1)
            locEn = λ * OsHE(j) + (1 - λ) * OsHB(j)
            A = MPO(locEn, sites)
            ComputedLocEn = real(inner(ψ', A, ψ)) - real(inner(ψ₀', A, ψ₀))
            En[n,j-1] = ComputedLocEn
            EnLog[n,j-1] = log10(abs(ComputedLocEn))
            Linkdims[n,j-1] = linkdim(ψ,j)
            Maxlinkdim[n] = maxlinkdim(ψ)
        end
        end

        # Print information
        println("Step ", n, " of ", N, " completed.")
        println("Maximum link dimension of MPS: ", maxlinkdim(ψ))
        println("Relative Energy: ", real(inner(ψ', H, ψ) - inner(ψ₀', H, ψ₀)))
        println("Estimated time remained: ", timeElapsed * (N - n), " seconds.")
    end

    # heatmap(positions, times, En)
    # savefig("simulation_output/En.png")

    # heatmap(positions, times, EnLog)
    # savefig("simulation_output/EnLog.png")

    # heatmap(positions, times, Linkdims)
    # savefig("simulation_output/Linkdims.png")

    titleString = "L = $L, λ = $λ, jⁱᵐᵖ = $jⁱᵐᵖ, jᵂᴾ = $jᵂᴾ, kᵂᴾ = $kᵂᴾ,
                    σᵂᴾ = $σᵂᴾ, ϵᵂᴾ = $ϵᵂᴾ, Nᴰᴹᴿᴳ = $Nᴰᴹᴿᴳ, χᴰᴹᴿᴳ ≤ $χᴰᴹᴿᴳ,
                    ϵᴰᴹᴿᴳ = $ϵᴰᴹᴿᴳ, χᵀᴰⱽᴾ ≤ $χᵀᴰⱽᴾ, ϵᵀᴰⱽᴾ = $ϵᵀᴰⱽᴾ, dt = $dt, Δt = $Δt"

    l = @layout [grid(2,2); b{0.2h}]

    title = plot(title = titleString, grid = false, showaxis = false)
    p1 = heatmap(positions, times, En, title = "Energy density")
    p2 = heatmap(positions, times, EnLog, title = "Log10 of Energy dens.")
    p3 = heatmap(positions, times, Linkdims, title = "Link dims")
    p4 = plot(times, Maxlinkdim, title = "Max link dims")
    plot(p1, p2, p3, p4, title, layout = l, dpi = 300, size=(1000, 1000))
    savefig("simulation_output/All.png")
end