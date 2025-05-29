function dispersion_relation(H::Hermitian{Matrix}, T::Matrix, L::Integer; kwargs...)
    # TODO : implement the dispersion relation for the dense case
end

function dispersion_relation(H::Hermitian{SparseMatrixCSC}, T::SparseMatrixCSC, L::Integer; nlevels::Integer = 2, kwargs...)
    dimH = size(H)[1]
    if dimH != size(T)[1]
        @error "Size of H and T must be tha same."
    end
    Idm = operator_identity(SparseMatrixCSC, size(H)[1])
    if HType == SparseMatrixCSC || TOperatorType == Matrix
        @warn "H sparse and T dense (non sense). Consider using T sparse, if possible."
        @warn "Converting H to a dense matrix."
        H = Matrix(H)
    end
    print("Computing the dispersion relation...")
    println("")
    nonhermiticity = norm(H * T - T * H)
    print("|HT-TH| = ", nonhermiticity)
    println("")
    # We compute the highest and lowest energy and we compute the range
    print("Computing the highest and lowest energy levels...")
    Emin, Emax = 0.0, 0.0 # Def: the lowest and highest energy levels
    if HType <: Hermitian
        eng0, _ = eigen(H)
        Emax, Emin = maximum(eng0), minimum(eng0)
    elseif HType <: SparseMatrixCSC
        val, _, _ = eigsolve(H, size(H, 1), 1, :SR, eltype(H); ishermitian = true)
        Emin = real(val[1])
        val, _, _ = eigsolve(H, size(H, 1), 1, :LR, eltype(H); ishermitian = true)
        Emax = real(val[1])
    end
    print("\r\u001b[2K")
    # We define the energy range and the energy shift
    ΔE = Emax - Emin # Def: the energy range.
    div_fraction = 10 # Def: the energy range divisor (see later)
    λ = Emax + ΔE / div_fraction # Def: the energy shift
    # Def: the projector to the momentum k (generating function)
    function P(k)
        proj = Idm
        for _ in 1:L-1
            proj = exp(-im * k) * T * proj + Idm
        end
        return proj
    end
    # We project the Hamiltonian to the eigenspace of every momentum k and we
    # compute the dispersion relation.
    print("Projecting and diagonalizing the Hamiltonian for each momentum k...")
    println("")
    disprelvec = Vector{BlochState}()
    for kl in 0:(L-1) # We iterate over all the possible momenta
        # Def: f = k / π (which is assumed to be a fraction)
        f = (2 * kl) // L
        f = f > 1 ? f - 2 : f
        k = π * f
        print("Step $(kl+1) out $L. Momentum: $f π ") # We print the momentum (keep it for debug)
        # Projecting
        Hk = deepcopy(H) # Def: the Hamiltonian to project to the momentum k   
        Hk -= λ * I # We shift down
        # Note for the following step: [H,T]=0 => P(k) can be applied only to the left
        Hk = (P(k) * Hk) / L # We project the Hamiltonian
        if HType <: Hermitian
            Hk = Hermitian(Hk)
        end
        Hk += λ * I # We shift up back
        # We compute the eigenvalues and eigenvectors of the projected Hamiltonian
        en, st = [], []
        if HType <: Hermitian
            en, st = eigen(Hk)
            st = collect(eachcol(st))
        elseif HType <: SparseMatrixCSC
            en, st, _ = eigsolve(Hk, size(H, 1), nlevels, :SR, ComplexF64; ishermitian = true)
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
    println("")
    return disprelvec
end

# MPO case => method is Tensor Networks
function dispersion_relation(H0::MPO; nlevels::Int = 3, kwargs...)
    error("Method disprel with Pure Tensor Newtork approach not implemented yet :(")
end

export dispersion_relation