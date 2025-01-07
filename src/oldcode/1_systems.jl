include("0_basics.jl")

# IMPORTANT PHYLOSLPHY: all the getters of the library (like hamiltonian) does not compute
# the value if it is not already computed. This is because the user may want to compute
# the value with a specific method. If the value is not computed, the function will return
# an error. The function dedicated to compute the value will be called with the same name
# but with the suffix `_compute` (like hamiltonian_compute), or there are dedicated functions
# to compute the value (like groundstate_energy).

# SYSTEMS

""" Structure for a 1-dimensional quantum many-body system.
    It is a composed of nested dictionaries.

    Data:
    - `parameters`: dictionary with the numerical parameters of the system.
        - `L`: number of sites.
        - `d`: local dimension.
        - `N`: dimension of the Hilbert space.
    - `operators`: dictionary with the operators of the system.
        - `H`: Hamiltonian operator (can be encoded as a function).
        - `T`: Translation operator (can be encoded as a function).
        - `Hj`: the local Hamiltonian operator at a certain site.
        - ... (define also the local operators).
    - `properties`: dictionary with the properties of the system.
        - `gap`: the energy gap of the system.
        - `disprel`: the dispersion relation of the system.
        - `translationinv`: if the system is translation invariant (true or false).
    """
mutable struct QuantumSystem
    data::Dict{Any, Any}
end

function QuantumSystem()
    S = QuantumSystem(Dict(
        :parameters => Dict(),
        :operators => Dict(),
        :states => Dict(),
        :properties => Dict()
    ))
    return S
end

""" Get the number of sites of a quantum system.
    """
function systemsize(S::QuantumSystem; kwargs...)::Int
    S.data[:parameters][:L]
end

""" Get the local dimension of a quantum system.
    """
function localdim(S::QuantumSystem; kwargs...)::Int
    S.data[:operators][:d]
end

""" Get the dimension of the Hilbert space of a quantum system.
    """
function hamiltonian(S::QuantumSystem; kwargs...)::QuantumOperator
    S.data[:parameters][:H]
end

""" Get the translation operator of a quantum system.
    """
function translation_operator(S::QuantumSystem; kwargs...)::QuantumOperator
    S.data[:parameters][:T]
end

# OPERATORS

""" The structure for a quantum operator.
    It is a dictionary with the data of the operator.
    
    Data:
    - `parentsystem`: the Quantum System(s) which it belongs to.
    - `representation`: a dictionary with structure {Type => Any}.
    - `eigenstates`: a list of QuantumStates.
    - `commuteswith`: a list of operators with which it commutes.
    - `ishermitian`: true or false.
    - `isunitary`: true or false.
    - `isidempotent`: true or false.
    - `supportsize`: the size of the support of the operator. 1 for local operators.
    """
mutable struct QuantumOperator
    data::Dict{Any, Any}
end

function get_parentsystem(O::QuantumOperator)::QuantumSystem
    O.data[:parentsystem]
end

function set_parentsystem!(S::QuantumSystem, O::QuantumOperator)::Nothing
    O.data[:parentsystem] = S
end

function get_representation(O::QuantumOperator; type::Type)::Any
    O.data[:representation][type]
end

function push_representation!(O::QuantumOperator, type::Type, value::Any)::Nothing
    O.data[:representation][type] = value
end

function get_spectrum(O::QuantumOperator)::Vector{Number}
    # We select all the values of the eigensystem
    sort(Set(values(O.data[:eigensystem])))
end

""" Compute the groundstate energy of a quantum system.

    # Arguments
    - `S::QuantumSystem`: the quantum system.

    # Keyword arguments
    - `which`: It is a symbol or a list of symbols. It specifies which eigenvalues are wanted.
        The possible values are:
        `:all` (default) tries to compute all the eigenvalues,
        `:min` computes only the minimum energy,
        `:max` computes only the maximum energy,
        `:maxmin` computes the maximum and minimum energy,
        `:lowestn` computes the first `n` eigenvalues (with minimum energy, in that case
        `n` must be provided as a keyword argument),
        `:highestn` computes the last `n` eigenvalues (with maximum energy, in that case
        `n` must be provided as a keyword argument).
    - `nlevels::Int`: The number of levels to compute. Default is 1.
    - `method::Symbol`: The method to compute the groundstate.
        Don't provide it or set `:auto` to use the default method.
    - `storestates::Bool`: If also the eigenstates should be stored. Default is `true`.
    """
function compute_spectrum(O::QuantumOperator; kwargs...)
    R = representation(O)
    if kwargs[:which] == :all && typeof(R) == Hermitian && kwargs[:storestates] == true
        λ, v = eigen(R)
        O.data[:maxeigenval] = λ[end]
        O.data[:mineigenval] = λ[1]
        O.data[:spectrum] = λ
        for vec in v
            state = QuantumState()
            if haskey(O.data, :parentsystem)
                state.data[:system] = O.data[:parentsystem]
            end
            state.data[:wavefunction] = vec
            push!(O.data[:eigenstates], state)
        end
    end
    if kwargs[:which] == :maxmin && typeof(R) == Hermitian && kwargs[:storestates] == false
        λ, _ = eigen(R)
        O.data[:maxeigenval] = λ[end]
        O.data[:mineigenval] = λ[1]
        O.data[:spectrum] = [λ[1], λ[end]]
    else
        error("Arguments not supported.")
    end
end

function eigenstates(O::QuantumOperator)
    keys(O.data[:eigensystem])
end

function eigenvalues(O::QuantumOperator)
    values(O.data[:eigensystem])
end

function ishermitian(O::QuantumOperator)::Bool
    O.data[:ishermitian]
end

function isunitary(O::QuantumOperator)::Bool
    O.data[:isunitary]
end

function isidempotent(O::QuantumOperator)::Bool
    O.data[:isidempotent]
end

function supportsize(O::QuantumOperator)::Int
    O.data[:supportsize]
end

function islocal(O::QuantumOperator)::Bool
    supportsize(O) == 1
end

function maxeigenval(O::QuantumOperator)
    O.data[:maxeigenval]
end

function mineigenval(O::QuantumOperator)
    O.data[:mineigenval]
end

# STATES

"""
    Structure for a quantum state.

    Naming conventions:
    - `name`: an optional name associated to the state.
    - `parentsystem`: the system associated to the state.
    - `wavefunction`: the wavefunction representation of the state (any form).
    - `quantumnumber`: a dictionary of structure {Operator => Eigenvalue}. 
    - `properties`: a dictionary with the properties of the state.
        - `E`: the average energy of the state.
        - `k`: the average momentum of the state.
    """
mutable struct QuantumState
    data::Dict{Any, Any}
end

function energy(state::QuantumState)
    state.data[:quantumnumbers][:E]
end

function momentum(state::QuantumState)
    state.data[:quantumnumbers][:k]
end

function average_energy(state::QuantumState)
    state.data[:properties][:E]
end

function average_momentum(state::QuantumState)
    state.data[:properties][:k]
end

# MODELS

""" Generates a quantum Ising model chain system.
    The Hamiltonian operator of the system reads:

    ```math
    H = -Jᶻ ∑ⱼ σᶻⱼσᶻⱼ₊₁ + hˣ ∑ⱼ σˣⱼ.
    ```

    # Arguments
    - `L::Int`: The size of the chain.

    # Optional arguments
    - `type::Type`: The type of the Hamiltonian. Default is `Matrix`.

    # Keyword arguments
    - `Jx::Number`: The coupling constant between spins in the x-direction.
    - `Jz::Number`: The coupling constant between spins in the z-direction.
    - `hx::Number`: The external field in the x-direction.
    - `hz::Number`: The external field in the z-direction.
    - `BC::Symbol`: The boundary conditions. Set `:pbc` for periodic boundary conditions.
    """
function isingmodel(L::Int; type::Type = Matrix, kwargs...)
    # We initialize the system
    S = QuantumSystem(Dict(), Dict())

    # We give a name to the system
    S.parameters[:name] = "Ising Model"

    # We store the number of sites
    S.parameters[:L] = L

    # We store the spin
    spin = 1//2
    S.parameters[:spin] = spin

    # We store the local dimension
    d = Int(2 * spin + 1)
    S.parameters[:d] = d

    # We store the dimension of the Hilbert space
    N = d^L
    S.parameters[:N] = N

    if type == Matrix

        # We define the local operators
        σx = Matrix([0 1; 1 0])
        σz = Matrix([1 0; 0 -1])

        # We define the translation operator
        T = generate_translation_operator_matrix(2, L)
        S.operators[:T] = T

        # We define the product of local operators as a shortcut
        𝒪(args...) = product_locals(L, 2, args...)

        if haskey(kwargs, :BC)
            BC = kwargs[:BC]
            S.parameters[:BC] = BC
        end

        # We define the Hamiltonian
        H = zeros(N, N)
        S.operators[:H] = H

        # We add a x-x interaction term if it is provided
        if haskey(kwargs, :Jx)
            Jx = kwargs[:Jx]
            S.parameters[:Jx] = Jx
            S.operators[:H] += sum(𝒪((-Jx * σx,j),(σx,j+1)) for j in 1:L-1)
            if S.parameters[:BC] == :pbc
                S.operators[:H] += 𝒪((-Jx * σx,L),(σx,1))
            end
        end

        # We add a z-z interaction term if it is provided
        if haskey(kwargs, :Jz)
            Jz = kwargs[:Jz]
            S.parameters[:Jz] = Jz
            S.operators[:H] += sum(𝒪((-Jz * σz,j),(σz,j+1)) for j in 1:L-1)
            if S.parameters[:BC] == :pbc
                S.operators[:H] += 𝒪((-Jz * σz,L),(σz,1))
            end
        end

        # We add a z-field if it is provided
        if haskey(kwargs, :hz)
            hz = kwargs[:hz]
            S.parameters[:hz] = hz
            S.operators[:H] += sum(𝒪((hz * σz, j)) for j in 1:L)
        end

        # We add a x-field if it is provided
        if haskey(kwargs, :hx)
            hx = kwargs[:hx]
            S.parameters[:hx] = hx
            S.operators[:H] += sum(𝒪((hx * σx, j)) for j in 1:L)
        end
    else
        error("The type is not supported.")
    end

    return S
end

""" Generates a quantum Clock model chain system.
    The Hamiltonian operator of the system reads:

    ```math
    H = - ∑ⱼ (J σⱼ⁺σⱼ₊₁ + H.c.) - ∑ⱼ (f τⱼ + H.c.) - ∑ⱼ (h σⱼ + H.c.) ,
    ```

    where the couplings can be complex (chiral).

    # Arguments
    - `L::Int`: The size of the chain.

    # Optional arguments
    - `n::Int`: The number of local states. Default is 3.
    - `type::Type`: The type of the Hamiltonian. Default is `Matrix`.

    # Keyword arguments
    - `J::Number`: The coupling constant between spins. Includes chirality.
    - `f::Number`: The transverse field (τ) strength. Includes chirality.
    - `h::Number`: The longitudinal field (σ) strength. Includes chirality.
    - `BC::Symbol`: The boundary conditions. Set `:pbc` for periodic boundary conditions.
    """
function clockmodel(L::Int; n::Int = 3, type::Type = Matrix, kwargs...)

    S = QuantumSystem(Dict(), Dict())

    # We give a name to the system
    S.parameters[:name] = "Quantum Clock Model"

    # We define the main parameters
    S.parameters[:L] = L # The number of sites

    # We check that n is a positive integer >= 2
    if !isinteger(n) || n < 2
        error("The clock number n must be a positive integer >= 2.")
    end

    # We store the clock number
    S.parameters[:n] = n

    # We store the local dimension
    d = n
    S.parameters[:d] = n

    # We store the dimension of the Hilbert space
    N = d^L
    S.parameters[:N] = N

    # We define the unit of phase
    ω = exp((2π * im) / n)

    if type == Matrix

        # We define the local operator σ
        σ = zeros(ComplexF64, n, n)
        for k in 1:n
            σ[k,k] = ω^(k-1)
        end

        # We define the local operator τ
        τ = zeros(ComplexF64, n, n)
        for k in 2:n
            σ[k,k-1] = 1
        end
        τ[1,n] = 1

        # We define the translation operator
        T = generate_translation_operator_matrix(d, L)
        S.operators[:T] = T

        # We define the product of local operators as a shortcut
        𝒪(args...) = product_locals(L, d, args...)

        # We store the boundary conditions
        if haskey(kwargs, :BC)
            BC = kwargs[:BC]
            S.parameters[:BC] = BC
        end

        # We add a σ-σ interaction term if it is provided
        if haskey(kwargs, :J)
            J = kwargs[:J]
            S.parameters[:J] = J
            HJ += sum(𝒪((-J * σ',j),(σ,j+1)) for j in 1:L-1)
            if S.parameters[:BC] == :pbc
                HJ += 𝒪((-J * σ',L),(σ,1))
            end
            S.operators[:H] += HJ + HJ'
        end

        # We add a σ-field term if it is provided
        if haskey(kwargs, :h)
            h = kwargs[:h]
            S.parameters[:h] = h
            Hh += sum(𝒪((-h * σ, j)) for j in 1:L)
            S.operators[:H] += Hh + Hh'
        end

        # We add a σ-field term if it is provided
        if haskey(kwargs, :f)
            f = kwargs[:f]
            S.parameters[:f] = f
            Hf += sum(𝒪((-f * τ, j)) for j in 1:L)
            S.operators[:H] += Hf + Hf'
        end

        # We add a longitudinal field if it is provided
        if haskey(kwargs, :hσ)
            hσ = kwargs[:hσ]
            H += sum(𝒪((-hσ * σ, j)) for j in 1:L)
        end
    else
        error("The type is not supported.")
    end

    return S
end

""" Generates a quantum Bose-Hubbard model chain system.
    The Hamiltonian operator of the system reads:

    ```math
    H = -t ∑ (bⱼ⁺bⱼ₊₁ + H.c.) + U/2 ∑ⱼ nⱼ(nⱼ-1) - μ ∑ⱼ nⱼ + V ∑ⱼ nⱼnⱼ₊₁.
    ```

    and it is truncated to the local dimension d.
    
    # Arguments
    - `L::Int`: The Hamiltonian of the system.

    # Optional arguments
    - `d::Int`: The local dimension, which is the truncation of the local space.
                Default is 2.
    - `type::Type`: The type of the Hamiltonian. Default is `Matrix`.

    # Keywords Arguments
    - `t::Real`: The hopping parameter.
    - `U::Real`: The on-site interaction parameter.
    - `μ::Real`: The chemical potential.
    - `V::Real`: The nearest-neighbor interaction parameter.
    - `BC::Symbol`: The boundary conditions.
                    Default is `:obc` for open boundary conditions.
    """
function bosehubbard(L::Int; d::Int = 2, type::Type = Matrix, kwargs...)

    S = QuantumSystem(Dict(), Dict())

    # We give a name to the system
    S.parameters[:name] = "Bose-Hubbard model"

    # We define the main parameters
    S.parameters[:L] = L # The number of sites

    # We check that n is a positive integer >= 2
    if !isinteger(d) || d < 2
        error("The local dimension d must be a positive integer >= 2.")
    end

    # We store the local dimension
    S.parameters[:d] = d

    # We store the dimension of the Hilbert space
    N = d^L
    S.parameters[:N] = N

    if type == Matrix

        # We define the local operator n
        n = zeros(ComplexF64, d, d)
        for k in 1:d
            σ[k,k] = k-1
        end

        # We define the local operator b
        b = zeros(ComplexF64, d, d)
        for k in 2:n
            σ[k,k-1] = √(k-1)
        end

        # We define the translation operator
        T = generate_translation_operator_matrix(d, L)
        S.operators[:T] = T

        # We define the product of local operators as a shortcut
        𝒪(args...) = product_locals(L, d, args...)

        # We store the boundary conditions
        if haskey(kwargs, :BC)
            BC = kwargs[:BC]
            S.parameters[:BC] = BC
        end

        # -t ∑ (bⱼ⁺bⱼ₊₁ + H.c.)
        # We add the hopping term if it is provided
        if haskey(kwargs, :t)
            t = kwargs[:t]
            S.parameters[:t] = t
            Ht += sum(𝒪((-t * b',j),(b,j+1)) for j in 1:L-1)
            if S.parameters[:BC] == :pbc
                Ht += 𝒪((-t * b',L),(b,1))
            end
            S.operators[:H] += Ht + Ht'
        end

        # + U/2 ∑ⱼ nⱼ(nⱼ-1)
        # We add the on-site interaction term if it is provided
        if haskey(kwargs, :U)
            U = kwargs[:U]
            S.parameters[:U] = U
            HU += sum(𝒪((U/2 * n * (n - I), j)) for j in 1:L)
            S.operators[:H] += HU
        end

        # - μ ∑ⱼ nⱼ
        # We add the chemical potential term if it is provided
        if haskey(kwargs, :μ)
            μ = kwargs[:μ]
            S.parameters[:μ] = μ
            Hμ += sum(𝒪((-μ * n, j)) for j in 1:L)
            S.operators[:H] += Hμ
        end

        # + V ∑ⱼ nⱼnⱼ₊₁
        # We add the nearest-neighbor interaction term if it is provided
        if haskey(kwargs, :V)
            V = kwargs[:V]
            S.parameters[:V] = V
            HV += sum(𝒪((V * n, j),(n, j+1)) for j in 1:L-1)
            if S.parameters[:BC] == :pbc
                HV += 𝒪((V * n, L),(n, 1))
            end
            S.operators[:H] += HV
        end
    else
        error("The type is not supported.")
    end

    return S
end

""" Computes the excitations of a quantum System.

    ## Arguments
    - `S::QuantumSystem` is the quantum system.

    ## Optional arguments
    - `k::Real` is the momentum of the excitation. If not provided, all the
       momenta are computed.
    - `nlevels::Int` is the number of excitations (the lowest ones).
    - `method::Symbol` is the method to compute the excitations. One can choose between,
      `:exact` (Exact diagonalization), `:krylov` (Lanczos / Arnoldi),
      `:tns` (Tensor Networks). If not provided, the method is chosen automatically,
       which is same of giving `:auto`.
    """

function excitations(
        S::QuantumSystem;
        k::Union{Real, Nothing} = nothing,
        nlevels::Int = 1,
        method::Symbol = :auto,
        kwargs...
    )::Set{QuantumState}

    # We get the parameters of the system, obtaining error if not defined
    L = systemsize(S)
    H = deepcopy(hamiltonian(S))
    Hk = deepcopy(H)
    T = translation_operator(S)
    Emin = energy(groundstate(S))
    Emax = energy(higheststate(S))
    ΔE = abs(Emax - Emin)

    # We divide the energy range in M parts. The parameter λ is a shift in the Hamiltonian
    # such that all the eigenvalues are for sure negative. This allows to compute the
    # eigenspace of the Hamiltonian after projection of the momentum k. Hence, M is
    # a conventional parameter, and we set it to 10.
    M = 10
    λ = Emax + ΔE/M

    if method == :exact

        # We assert that the Hamiltonian is translational invariant
        if norm(H * T - T * H) > 1e-10
            error("The Hamiltonian is not translational invariant.")
        end

        # We project the Hamiltonian to the eigenspace of momentum k.
        if (H isa AbstractMatrix) && (H isa AbstractMatrix)
            H = H - λ * I
            Hk = Hk - λ * I
            for _ in 1:(L-1)
                Hk = H + exp(-im * k) * T * Hk
            end
            Hk /= L
            Hk =  Hk + λ * I
        else
            error("Types of H and T not supported for the function.")
        end
        
        set::Set{QuantumState} = Set()

        energies, states = eigen(Hermitian(Matrix(Hk)))

        for i in 1:nlevels
            if i > length(energies) || energies[i] > Emax + ΔE/(2M)
                break
            end
            push!(set, QuantumState(S, states[:,i], Dict(:E => energies[i])))
        end

        # We add the momentum
        for state in set
            state.properties[:k] = k
        end

        return set

    elseif method == 

    # If S has not L, throw error
    if !haskey(S.parameters, :L)
        error("The number of sites L is not defined in the system.")
    else
        L = S.parameters[:L]
    end

    # If S has not H, throw error
    if !haskey(S.operators, :H)
        error("The Hamiltonian operator is not defined in the system.")
    else
        H0 = deepcopy(hamiltonian(S))
    end

    # If S has not T, throw error
    if !haskey(S.operators, :T)
        error("The translation operator is not defined in the system.")
    else
        T = S.operators[:T]
    end

    # We project the Hamiltonian to the eigenspace of momentum k.
    if (H0 isa AbstractMatrix) && (H0 isa AbstractMatrix)
        for _ in 1:(L-1)
            H0 += exp(im * k) * T * H0
        end
        H0 /= L
    else
        error("Types of H and T not supported for the function.")
    end
    
    set::Set{QuantumState} = Set()

    # Compute the eigenvalues of Hk, and we sleect the first n values and vectors.
    # We return the energies and the first n states.
    energies, states, info = eigsolve(H0; howmany = nlevels, which = :SR, T = ComplexF64)

    for i in 1:nlevels
        push!(set, QuantumState(S, states[:,i], Dict(:E => energies[i])))
    end

    # We add the momentum
    for state in set
        state.properties[:k] = k
    end

    return set
end

function disprel_exact(
        S::QuantumSystem, # The quantum system
        nlevels::Int, # The number of excitations (lowest ones)
    )::Set{QuantumState}

    # If S has not L, throw error
    if !haskey(S.parameters, :L)
        error("The number of sites L is not defined in the system.")
    else
        L = S.parameters[:L]
    end

    set::Set{QuantumState} = Set()

    for m in 0:L-1
        f = m / L
        if f > 0.5
            f = f - 1
        end
        k = 2π * f
        e = excitations_exact(S, k, nlevels)
        for state in e
            push!(set, state)
        end
    end

    return set
end

function plot_disprel(set::Set{QuantumState}; aspect_ratio = 1)
    k = [momentum(state) for state in set]
    E = [energy(state) for state in set]

    plot(k, E, seriestype = :scatter,
        legend = false,
        aspect_ratio = aspect_ratio,
        xlims = (-π, π),
        xticks = (-π:π/2:π, ["-π", "-π/2", "0", "π/2", "π"]),
        dpi = 300)
end

S = clockmodel(3, 7, 1.0 + 0im, hτ = 1.0 * exp(im * 0/3), hσ = 0.2 + 0im, BC = :pbc)
# S = isingmodel(7, 1; hx = -1, BC = :pbc)
setstates = disprel_exact(S, 100000)
plot_disprel(setstates, aspect_ratio = 1)