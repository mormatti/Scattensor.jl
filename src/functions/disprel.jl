# DEFINITION
# A generic vector of bloch states can be a whole dispersion relation, a band, etc.

# PLOTTING

function Plots.plot(disprelvec::Vector{T}
    ) where {T <: BlochState}

    E = [energy(state) for state in disprelvec]
    k = [momentum(state) for state in disprelvec]
    return Plots.plot(k, E, seriestype = :scatter, label = "Dispersion relation")
end

# PROPERTIES

# TODO: find a way to write these function in a more efficient 
# (complexity) and compact way

function getgroundstate(states::Vector{BlochStateType}
    ) where {BlochStateType <: BlochState}
    # We identify the position of the state which have the lowest energy
    i0 = 1
    for i in eachindex(states)
        if energy(states[i]) < energy(states[i0])
            i0 = i
        end
    end
    return states[i0]
end
export getgroundstate

function popgroundstate!(states::Vector{BlochStateType}
    ) where {BlochStateType <: BlochState}
    # We identify the position of the state which have the lowest energy
    i0 = 1
    for i in eachindex(states)
        if energy(states[i]) < energy(states[i0])
            i0 = i
        end
    end
    return popat!(states, i0)
end
export popgroundstate!

function getstatesabove(states::Vector{BlochStateType}, energy::RealType
    ) where {BlochStateType <: BlochState, RealType <: Real}
    return filter(state -> energy(state) <= energy, states)
end
export getstatesabove

function getstatesbelow(states::Vector{BlochStateType}, energy::RealType
    ) where {BlochStateType <: BlochState, RealType <: Real}
    return filter(state -> energy(state) >= energy, states)
end
export getstatesbelow

function selectfirstband(states::Vector{BlochStateType}
    ) where {BlochStateType <: BlochState}

    if isempty(states)
        return Vector{BlochStateType}()
    else
        band = [getgroundstate(states)]
        pop!(band)
        while !isempty(states)
            ground = getgroundstate(states)
            k = momentum(ground)
            # We remove all the states which have momentum k
            states = filter(state -> momentum(state) != k, states)
            push!(band, ground)
        end
        return band
    end
end
export selectfirstband

function generateidentity(type::Type, n::IntType) where {IntType <: Int}
    if type <: Hermitian
        return Hermitian(Matrix{Float64}(I, n, n))
    elseif type <: Matrix
        return Matrix{Float64}(I, n, n)
    elseif type <: SparseMatrixCSC
        return sparse(I, n, n)
    else
        error("Type $(type) not supported in function generateidentity.")
    end
end
export generateidentity

# CONSTRUCTORS

function disprel(H0::LocalOperator{MatrixType}, L::Int;
    nlevels::Int = 2, 
    kwargs...
    ) where {MatrixType <: Union{Hermitian, SparseMatrixCSC}}
    
    # We construct the Hamiltonian and the translation operator
    if !hasuniformlocaldims(H0)
        error("The local dimensions of the sites must be the same.")
    end
    d = H0.localdims[1]
    L0 = length(H0) # Def: the number of sites of the local Hamiltonian
    Idm = generateidentity(MatrixType, d^L) # Def: the identity matrix for the local space
    Idmext = generateidentity(MatrixType, d^(L-L0)) # Def: the identity matrix for the complementary space
    H0ext = kron(H0.repr, Idmext)
    T = generate_translation_operator_matrix(d, L) # We obtain the translation operator T
    # We compute H, the Hamiltonian of the chain of length L
    H = deepcopy(H0ext)
    for _ in 1:L-1
        H0ext = T * H0ext * T'
        H += H0ext
    end
    if MatrixType <: Hermitian
        H = Hermitian(H)
    end

    # We compute the highest and lowest energy and we compute the range
    @info "Computing the highest and lowest energy..."
    Emin = 0.0 # Def: the lowest energy
    Emax = 0.0 # Def: the highest energy
    if MatrixType <: Hermitian
        eng0, _ = eigen(H)
        Emax, Emin = maximum(eng0), minimum(eng0)
    elseif MatrixType <: SparseMatrixCSC
        val, _, _ = eigsolve(H, size(H, 1), 1, :SR, eltype(H); ishermitian = true)
        Emin = real(val[1])
        val, _, _ = eigsolve(H, size(H, 1), 1, :LR, eltype(H); ishermitian = true)
        Emax = real(val[1])
    end
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
    @info "Projecting and diagonalizing the Hamiltonian for each momentum k..."
    disprelvec = Vector{BlochState}()
    print("Momentum: ")
    for kl in 0:(L-1) # We iterate over all the possible momenta

        print("$(kl)π ") # We print the momentum (keep it for debug)

        # Def: f = k / π (which is assumed to be a fraction)
        f = (2 * kl) // L
        f = f > 1 ? f - 2 : f
        k = π * f

        # Projecting
        Hk = deepcopy(H) # Def: the Hamiltonian to project to the momentum k   
        Hk -= λ * I # We shift down
        # Note for the following step: [H,T]=0 => P(k) can be applied only to the left
        Hk = (P(k) * Hk) / L # We project the Hamiltonian
        if MatrixType <: Hermitian
            Hk = Hermitian(Hk)
        end
        Hk += λ * I # We shift up back

        # We compute the eigenvalues and eigenvectors of the projected Hamiltonian
        en, st = [], []
        if MatrixType <: Hermitian
            en, st = eigen(Hk)
            st = collect(eachcol(st))
        elseif MatrixType <: SparseMatrixCSC
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
    end

    return disprelvec
end

# MPO case => method is Tensor Networks
function disprel(H0::LocalOperator{MPO};
                    nlevels::Int = 3,
                    kwargs...)
    
    error("Method disprel with Pure Tensor Newtork approach not implemented yet :(")
end