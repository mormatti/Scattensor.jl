# Basic utils

"""A shortcut binary notation for the periodic modulus."""
function ↻(n::Integer, m::Integer)::Integer
    n > 0 ? (n-1)%m + 1 : m + n%m
end

"""A shortcut binary notation for the periodic modulus centered in zero"""
function ZZ(n::Integer, p::Integer)::Integer
    s = ((n+1)/2)
    p = p + s
    return (p > 0 ? (p-1)%n + 1 : n + p%n) - s
end

"""A function to print colored text in the standard output."""
function printcolored(r, g, b, text)
    print("\e[1m\e[38;2;$r;$g;$b;249m", text)
end

# Matrices

""" Generates the whole matrix corresponding to the action of a local operators product.

    ## Inputs
    - `L` is the number of sites of the chain;
    - `d` is the local dimension;
    - `args` is a list of pairs (𝐚, j) where 𝐚 is the local operator written in the local 
    space (small matrix) and j is the position of the local operator 𝐚.
    """    
function product_locals(L::Int, d::Int, args::Vararg{Tuple{AbstractMatrix, Int}})
    N = d^L
    
    # Ensure all matrices match the local space dimension
    for (op, _) in args
        @assert size(op) == (d, d) "Each local operator must have dimensions $d x $d."
    end
    
    # Normalize positions modulo L and store them
    args = [(op, pos ↻ L) for (op, pos) in args]
    
    # Identity matrix for the local space
    𝟙 = SparseMatrixCSC(Matrix{Int}(I, d, d))
    
    # Build operators list, filling with identity when no operators are specified
    𝒱 = [reduce(*, [op for (op, pos) in args if pos == i]; init = 𝟙) for i in 1:L]
    
    # Tensor product of the L site operators
    𝐀 = reduce(kron, 𝒱)
    @assert size(𝐀) == (N, N) "The resulting matrix size is incorrect."
    
    return 𝐀
end

""" Generates the translation operator for a chain of `L` sites with local dimension `d`.
    
    ## Assumptions
     - The system is assumed to be uniform, i.e. the local dimension is the same for all sites.
     - The system, in order to perform a translation, must be in periodic boundary conditions.

    ## Inputs
    - `L` is the number of sites of the chain.
    - `d` is the local dimension.

    ## Outputs
    - The translation operator `T` in matrix form.
    """
function generate_translation_operator_matrix(d::Integer, L::Integer)::SparseMatrixCSC{Int64, Int64}  
    N = d^L
    T = spzeros(Float64, N, N)

    Lst = []
    c = 0
    for _ in 1:d
        lst = []
        for _ in 1:(N/d)
            c = c + 1
            push!(lst, c)
        end
        push!(Lst, lst)
    end

    for indL in eachindex(Lst)
        lst = Lst[indL]
        for ind in eachindex(lst)
            j = lst[ind]
            T[j, ((d*(j-1)+1)%N) + indL - 1] = 1
        end
    end

    return T
end

# Useful but old function...
""" Simultaneously diagonalize an Hermitian matrix 𝐇 and a unitary matrix 𝐔 such that they commute, i.e. [𝐇,𝐔] = 0.

    ## Assumptions
    - `𝐇` must be translationally invariant, i.e. [𝐇,𝐔] = 0.

    ## Inputs
    - `𝐇` is an Hermitian matrix, i.e. 𝐇 = 𝐇†;
    - `𝐔` is a unitary matrix, i.e. 𝐔𝐔† = 𝐔†𝐔 = 𝟙.

    ## Outputs
    - `u` is the `Vector` of the phases (angles) of the eigenvalues of 𝐔;
    - `h` is the `Vector` of real eigenstates associated to 𝐇;
    - `𝛙` is the `Vector` of (common) eigenvectors associated to 𝐇 and 𝐔.
    """
function simultaneous_diagonalization_HU(
    𝐇::Matrix{ComplexF64},
    𝐔::Matrix{ComplexF64};
    check_hermiticity::Bool = false,
    check_unitarity::Bool = false,
    check_translational_invariance::Bool = false
    )::Tuple{Vector{Float64}, Vector{Float64}, Vector{Vector{ComplexF64}}}

    @debug "Computing the bloch states from exact diagonalization..."

    # We assert the two matrices have the same size and are square
    @assert size(𝐇)[1] == size(𝐇)[2] == size(𝐔)[1] == size(𝐔)[2] "𝐇 and 𝐔 must be square and having the same size."

    N = size(𝐇)[1] # The dimension of the Hilbert space

    if check_hermiticity
        @assert ishermitian(𝐇) "𝐇 is not hermitian. |𝐇† - 𝐇|/|𝐇| = $(norm(𝐇' - 𝐇)/norm(𝐇))."
    end 

    if check_unitarity
        @assert 𝐔 * 𝐔' ≈ 𝐔' * 𝐔 ≈ I "𝐔 is not unitary. |𝐔†𝐔 - I| = $(norm(𝐔' * 𝐔 - I))."
    end

    if check_translational_invariance
        @assert 𝐇 * 𝐔 ≈ 𝐔 * 𝐇 "𝐇 is not translational invariant. |𝐇𝐔 - 𝐔𝐇|/|𝐇𝐔| = $(norm(𝐇𝐔 - 𝐔𝐇)/norm(𝐇𝐔))."
    end

    # We compute the groundstate energy E₀ and the matrix product 𝐇𝐓,
    # where 𝐇 is shifted by E0 - 1
    h₀ = eigen(𝐇, permute = true).values[1]
    𝐇′𝐔 = (𝐇 - h₀ * I + I) * 𝐔

    # We compute all the eigenvectors and eigenvalues of HT
    (𝜆, 𝛙) = eigen(𝐇′𝐔)
    Nᵥ = length(𝜆) # The number of eigenvalues

    # Extract the phases and moduli from the eigenvalues 𝜆
    (u, h) = (zeros(Nᵥ), zeros(Nᵥ))
    for i in eachindex(𝜆)
        (u[i], h[i]) = (angle(𝜆[i]), abs(𝜆[i]) + h₀ - 1)
    end

    𝛙 = [𝛙[:,i] for i in 1:N]
    return u, h, 𝛙
end

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