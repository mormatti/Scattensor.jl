"""
    dispersion_relation(H::AbstractMatrix, T::AbstractMatrix, L::Integer; nlevels::Integer=2, kwargs...) -> Vector{BlochState}

Compute a dispersion relation by projecting a translation-invariant Hamiltonian onto momentum sectors.

Given a Hamiltonian `H` and a (unitary) translation operator `T` satisfying `[H, T] = 0`, this
function constructs projectors onto momentum eigenspaces and diagonalizes the projected Hamiltonian
for each allowed momentum `k = 2π m/L`.

# Arguments
- `H`: Hamiltonian matrix (dense or sparse), assumed Hermitian.
- `T`: Translation operator matrix, assumed unitary and of the same size as `H`.
- `L`: Number of lattice sites (sets the allowed momenta).

# Keyword Arguments
- `nlevels::Integer=2`: Number of lowest-energy eigenstates to keep per momentum sector.
- `kwargs...`: Reserved for future options (currently unused).

# Returns
- A vector of [`BlochState`](@ref) objects. Each element stores:
  - the eigenvector (as `data`),
  - the corresponding energy,
  - `koverpi = k/π` as a `Rational`.

# Side effects
- Prints diagnostic checks for Hermiticity/unitarity/translation invariance and progress messages.

# Notes
- Only nonnegative reduced-zone momenta are currently included (`koverpi ≥ 0`), so you may want to
  reconstruct `-k` points by symmetry when plotting/processing.
"""
function dispersion_relation(H::AbstractMatrix, T::AbstractMatrix, L::Integer; nlevels::Integer = 2, kwargs...)
    DenseType = Matrix
    SparseType = SparseMatrixCSC

    # Square matrix check for H
    if size(H)[1] != size(H)[2]
        error("The Hamiltonian is not a square matrix.")
    end

    # Square matrix check for T
    if size(T)[1] != size(T)[2]
        error("The translation operator is not a square matrix.")
    end

    # Size comparison check for H,T
    if size(H)[1] != size(T)[1]
        error("The Hamiltonian must have the same size of the translation operator.")
    end
    dimH = size(H)[1]

    # Hermiticity check for H
    nonhermiticity = norm(H - H')
    println("Hermiticity test: |H - H'| = $nonhermiticity")

    # Definition of identity matrix
    Idm = operator_identity(SparseMatrixCSC, dimH)

    # First unitarity check for T
    nonunitarity1 = norm(T * T' - Idm)
    nonunitarity2 = norm(T * T' - T' * T)
    println("First unitarity test: |TT' - I| = $nonunitarity1")

    # Second unitarity check for T
    nonunitarity2 = norm(T * T' - T' * T)
    println("Second unitarity test: |[T,T']| = $nonunitarity2")

    # Translational invariance check for H,T
    nontranslationalinv = norm(H * T - T * H)
    println("Translational invariance test: |HT - TH| = $nontranslationalinv")

    # If H is sparse and T dense, H is converted to a dense matrix
    if typeof(H) <: SparseType && typeof(T) <: DenseType
        @warn "H sparse and T dense (non sense). Consider using T sparse, if possible. Converting H to a dense matrix."
        H = Matrix(H)
    end

    println("Computing the dispersion relation...")
    println("Computing the highest and lowest energy levels...")
    Emin, Emax = 0.0, 0.0 # Def: the lowest and highest energy levels
    if typeof(H) <: DenseType
        eng0, _ = eigen(Hermitian(H))
        Emax, Emin = maximum(eng0), minimum(eng0)
    elseif typeof(H) <: SparseType
        val, _, _ = eigsolve(H, size(H, 1), 1, :SR, eltype(H); ishermitian = true)
        Emin = real(val[1])
        val, _, _ = eigsolve(H, size(H, 1), 1, :LR, eltype(H); ishermitian = true)
        Emax = real(val[1])
    else
        error("Type $(typeof(H)) non supported in the dispersion relation algorithm.")
    end

    print(specialchar_cancel_line())

    # We define the energy range and the energy shift
    ΔE = Emax - Emin # Def: the energy range.
    div_fraction = 10 # Def: the energy range divisor (see later)
    λ = Emax + ΔE / div_fraction # Def: the energy shift

    # Def: the projector to the momentum k (generating function)
    function projctr(k)
        result = Idm
        for _ in 1:L-1
            result = exp(-im * k) * T * result + Idm
        end
        return result
    end
    
    # We project the Hamiltonian to the eigenspace of every momentum k and we
    # compute the dispersion relation.
    println("Projecting and diagonalizing the Hamiltonian for each momentum k...")
    disprelvec = Vector{BlochState}()
    for kl in 0:(L-1) # We iterate over all the possible momenta
        # Def: f = k / π (which is assumed to be a fraction)
        f = (2 * kl) // L
        f = f > 1 ? f - 2 : f
        if f < 0
            continue
        end
        k = π * f
        println("Step $(kl+1) / $L. Momentum: $f⋅π ") # We print the momentum (keep it for debug).
        # Projecting
        Hk = deepcopy(H) - λ * I # Def: the Hamiltonian to project to the momentum k, here we shift it down
        Hk = (projctr(k) * Hk) / L + λ * I # We project the Hamiltonian and we shift up back. Note that P can be applied only to the left due to [H,T] = 0.
        en, st = [], [] # We compute the eigenvalues and eigenvectors of the projected Hamiltonian
        if typeof(H) <: DenseType
            en, st = eigen(Hermitian(Hk))
            st = collect(eachcol(st))
        elseif typeof(H) <: SparseType
            en, st, _ = eigsolve(Hk, size(H, 1), nlevels, :SR, ComplexF64; ishermitian = true, krylovdim = max(nlevels * 2, nlevels + 20))
            en = real(en)
        end
        for i in eachindex(en)
            if en[i] > Emax + ΔE/(2 * div_fraction)
                break
            end 
            blochstate = BlochState(st[i], en[i], f)
            push!(disprelvec, blochstate)
        end
        print("\r\u001b[2K")
    end
    print("Diagonalization done.")
    return disprelvec
end

"""
    dispersion_relation(H::MPO; nlevels::Int=3, kwargs...) -> Vector{BlochState}

Experimental tensor-network version of [`dispersion_relation`](@ref) operating on an `MPO`.

This method attempts to obtain momentum-resolved low-energy states using a penalty construction and
repeated DMRG runs. It is currently **under construction** and may be slow or unreliable depending
on model and parameters.

# Arguments
- `H::MPO`: Hamiltonian as an MPO with uniform local dimension.

# Keyword Arguments
- `nlevels::Int=3`: Number of low-energy states to target per momentum.
- `kwargs...`: Reserved for future options.

# Returns
- A vector of [`BlochState`](@ref) objects whose `data` field is an `MPS`.

# Warnings
- This routine is marked as TODO in the source: it is not yet a pure/robust TN projection method.
"""
function dispersion_relation(H::MPO; nlevels::Int = 3, maxdim = [50,100,200,400,500], cutoff = 1e-10, tol = 1e-8, kwargs...)
    
    dmrgargs = (nsweeps = 1000, maxdim = maxdim, cutoff = [cutoff], observer = CustomObserver(tol))

    # Define parameters and operators
    sites = siteinds_main(H)
    L = length(H)
    d = get_uniform_localdim(H)
    T = operator_translation(MPO, d, L)
    replace_siteinds!(T, sites)
    Idm = operator_identity(MPO, d, L)
    replace_siteinds!(Idm, sites)

    # Hermiticity and translational invariance checks
    nonhermiticity = norm(H - adjoint_mpo(H))
    println("Hermiticity test: |H - H'| = $nonhermiticity")
    nontranslationalinv = norm(product(H, T) - product(T, H))
    println("Translational invariance test: |HT - TH| = $nontranslationalinv")

    println("Computing the dispersion relation...")
    println("Computing the highest and lowest energy levels...")
    Emin, Emax = 0.0, 0.0 # Def: the lowest and highest energy levels

    # Finding the lowest and highest with DMRG
    # TODO improve with precition (DMRG listener)
    Emin, _ = dmrg(H, random_mps(sites, linkdims = 10); dmrgargs...)
    Emax, _ = dmrg(-H, random_mps(sites, linkdims = 10); dmrgargs...)
    Emax = -Emax

    # We define the energy range and the energy shift
    ΔE = Emax - Emin # Def: the energy range.
    safepenalty = ΔE / (1 - cos(2π / L))
    println("Safe penalty = ", safepenalty)

    # Def: the cosine penalty term
    function cospenalty(k)
        P = -0.5 * exp(im * k) * T
        mpotoreturn = Idm + P + adjoint_mpo(P)
        return mpotoreturn
    end
    
    # We project the Hamiltonian to the eigenspace of every momentum k and we
    # compute the dispersion relation.
    println("Projecting and diagonalizing the Hamiltonian for each momentum k...")
    disprelvec = Vector{BlochState}()
    for kl in 0:(L-1) # We iterate over all the possible momenta
        # Def: f = k / π (which is assumed to be a fraction)
        f = (2 * kl) // L
        f = f > 1 ? f - 2 : f
        if f < 0
            continue
        end
        k = π * f
        println("Step $(kl+1) / $L. Momentum: $f⋅π ") # We print the momentum (keep it for debug).

        # We project the Hamiltonian using the penalty term
        Hk = H + safepenalty * cospenalty(k)

        en = []
        st = Vector{MPS}()
        for level in 1:nlevels
            println("level: $level")
            energ, state = dmrg(Hk, st, random_mps(sites, linkdims = 10), weight = 100.0; dmrgargs...)
            push!(st, state)
            push!(en, energ)
        end

        for i in eachindex(en)
            blochstate = BlochState(st[i], en[i], f)
            push!(disprelvec, blochstate)
        end
    end
    print("Computation done.")
    return disprelvec
end
