"""
S-matrix utilities (experimental).

The routines in this file implement (or prototype) real-space and momentum-space two-particle
S-matrix constructions on top of ITensor's `MPS`/`MPO` machinery.

**Status**: under construction. The algorithms and APIs in this file are not yet considered stable,
and parts of the implementation still rely on provisional helpers (e.g. `substitute_siteinds!`).
"""

# TODO: all of this code is not yet fully functional, it needs to be tested and debugged.
# TODO: this code still use substitute_siteinds and blind_product_inner... fix it.

"""
    smatrix_real_space(d, H, ψ0, w, t, N, χmax) -> Array{ComplexF64,4}

Prototype real-space S-matrix expansion for two-particle scattering (MPO/MPS workflow).

This routine builds two-particle asymptotic in/out states by applying a localized creation operator `w`
at different positions on top of a vacuum `ψ0`, then computes a truncated time-evolution/series expansion
to obtain the real-space S-matrix tensor `S[j1, j2, j1′, j2′]`.

# Arguments
- `d::Integer`: Local dimension.
- `H::MPO`: Hamiltonian MPO governing time evolution.
- `ψ0::MPS`: Reference vacuum/ground state.
- `w::MPO`: Localized creation operator (as an MPO).
- `t::Real`: Evolution time.
- `N::Integer`: Expansion order / number of terms (implementation-dependent).
- `χmax::Integer`: Maximum bond dimension used in `apply`.

# Returns
- `S::Array{ComplexF64,4}` with indices `(j1, j2, j1′, j2′)` over the allowed particle positions.

# Warnings
- This function is marked TODO in the source and is not yet fully validated.
"""
function smatrix_real_space(d::Integer, H::MPO, ψ0::MPS, w::MPO, t::Real, N::Integer, χmax::Integer)
    # TODO assert that the local dimension of the MPS is the same as the one of the MPO
    # TODO assert that ψ0 is the groundstate of the Hamiltonian H
    Lψ0 = length(ψ0)
    l = length(w)
    L = Lψ0 - l + 1 # We define the length of the particle positions
    # We construct the single-particle basis wj
    println("Constructing the creation operators for length L...")
    W = Array{MPO}(undef, L)
    for j in 1:L
        W[j] = insert_local(j - 1, deepcopy(w), Lψ0 - (j - 1) - l, d)
        substitute_siteinds!(W[j], ψ0) # FIXME substitution MPS -> MPO siteinds
    end
    println("Constructing the asimptotic out states...")
    ψout = Array{MPS}(undef, L, L)
    for j1 in 1:L
        for j2 in 1:L
            ψtoadd = deepcopy(ψ0)
            ψtoadd = apply(W[j1], ψtoadd; maxdim = χmax)
            ψtoadd = apply(W[j2], ψtoadd; maxdim = χmax)
            normalize!(ψtoadd)
            ψout[j1, j2] = ψtoadd
        end
    end

    println("Constructing the asimptotic in states...")
    ψin = deepcopy(ψout)
    Loffset = Int64(round(Lψ0/4))
    Lmin = 1 + Loffset
    Lmax = L - Loffset


    # We compute the matrix A
    println("Computing the smatrix element expansion in real space...")
    S = zeros(ComplexF64, L, L, L, L) # N+1 because we include the N = 0 term
    for n in 0:N
        Sn = zeros(ComplexF64, L, L, L, L)
        normn = 0
        for j1 in Lmin:Lmax
            for j2 in Lmin:Lmax
                if abs(j1 - j2) < l
                    continue
                end
                factor = 1
                if n > 0
                    factor = -im * t / n
                end
                ψin[j1, j2] = factor * apply(H, ψin[j1, j2]; maxdim = χmax)
                for j1′ in 1:L
                    for j2′ in 1:L
                        if abs(j1′ - j2′) < l
                            continue
                        end
                        term = inner(ψout[j1′, j2′], ψin[j1, j2])
                        S[j1, j2, j1′, j2′] += term
                        Sn[j1, j2, j1′, j2′] = term
                        normn += abs(term)^2
                    end
                end
            end
        end
        println("Norm of the $n-th term: $(sqrt(normn))")
    end
    return S
end

export smatrix_real_space

"""
    smatrix_element_momentum_space(S, k1, k2, k1′, k2′) -> Complex

Compute a momentum-space S-matrix element from a real-space tensor `S`.

Given `S[j1, j2, j1′, j2′]`, this function performs the discrete Fourier transform over all four indices
to obtain the scattering amplitude between incoming momenta `(k1, k2)` and outgoing momenta `(k1′, k2′)`.

# Arguments
- `S::Array`: Real-space S-matrix tensor with 4 indices.
- `k1, k2, k1′, k2′::Real`: Momenta (radians) used in the phase factors.

# Returns
- A complex amplitude.
"""
function smatrix_element_momentum_space(S::Array, k1::Real, k2::Real, k1′::Real, k2′::Real)
    result = 0
    L1 = size(S, 1)
    L2 = size(S, 2)
    L1′ = size(S, 3)
    L2′ = size(S, 4)
    normaliz = sqrt(L1*L2*L1′*L2′)
    for j1 in 1:L1
        for j2 in 1:L2
            for j1′ in 1:L1′
                for j2′ in 1:L2′
                    prefactor = exp(-im * k1 * j1) * exp(-im * k2 * j2)
                    prefactor *= exp(im * k1′* j1′) * exp(im * k2′ * j2′)
                    prefactor /= normaliz
                    result += prefactor * S[j1, j2, j1′, j2′]
                end
            end
        end
    end
    return result
end

"""
    smatrix_momentum_space(S, k1, k2, k1′, k2′) -> Complex

Alias for [`smatrix_element_momentum_space`](@ref).
"""
smatrix_momentum_space(S::Array, k1::Real, k2::Real, k1′::Real, k2′::Real) =
    smatrix_element_momentum_space(S, k1, k2, k1′, k2′)

export smatrix_momentum_space

"""
    smatrix_real_space_tdvp(d, H, ψ0, w, t, χmax) -> Array{ComplexF64,4}

Prototype real-space S-matrix computation using TDVP time evolution.

This routine builds two-particle states by applying `w` at different positions and uses TDVP to evolve
them for a time `t`, then overlaps with the corresponding out-states to form the S-matrix tensor.

# Warnings
- This function is experimental and under construction.
"""
function smatrix_real_space_tdvp(
    d::dType, # The local dimension of the system
    H::MPO, # The Hamiltonian of the system
    ψ0::ψ0Type, # The vacuum state
    w::wType, # The Wannier creation operator
    t::tType, # The time at which we want to compute the S-matrix element
    χmax::Int # The maximum bond dimension
    ) where {
        dType <: Integer,
        ψ0Type <: MPS,
        wType <: MPO,
        tType <: Real
        }

    # TODO assert that the local dimension of the MPS is the same as the one of the MPO
    # TODO assert that ψ0 is the groundstate of the Hamiltonian H

    Lψ0 = length(ψ0)
    l = length(w)
    L = Lψ0 - l + 1 # We define the length of the particle positions

    # We construct the single-particle basis wj
    println("Constructing the creation operators for length L...")
    W = Array{MPO}(undef, L)
    for j in 1:L
        W[j] = insert_local(j - 1, deepcopy(w), Lψ0 - (j - 1) - l, d)
        substitute_siteinds!(W[j], ψ0)
    end

    println("Constructing the asimptotic out states...")
    ψout = Array{MPS}(undef, L, L)
    for j1 in 1:L
        for j2 in 1:L
            ψtoadd = deepcopy(ψ0)
            ψtoadd = apply(W[j1], ψtoadd; maxdim = χmax)
            ψtoadd = apply(W[j2], ψtoadd; maxdim = χmax)
            normalize!(ψtoadd)
            ψout[j1, j2] = ψtoadd
        end
    end

    println("Constructing the asimptotic in states...")
    ψin = deepcopy(ψout)
    Loffset = Int64(round(Lψ0/4))
    Lmin = 1 + Loffset
    Lmax = L - Loffset

    # We compute the matrix A
    println("Computing the smatrix element expansion in real space...")
    S = zeros(ComplexF64, L, L, L, L) # N+1 because we include the N = 0 term
    for j1 in Lmin:Lmax
        for j2 in Lmin:Lmax
            println("step $j1, $j2")
            if abs(j1 - j2) < l
                continue
            end
            dt = t / 10
            for n in 1:10
                ψin[j1, j2] = tdvp(H, -im * dt, ψin[j1, j2]; maxdim = χmax)
                println("iteration $n, maxlinkdim = $(maxlinkdim(ψin[j1, j2]))")
            end
            for j1′ in 1:L
                for j2′ in 1:L
                    if abs(j1′ - j2′) < l
                        continue
                    end
                    term = inner(ψout[j1′, j2′], ψin[j1, j2])
                    S[j1, j2, j1′, j2′] += term
                end
            end
        end
    end
    return S
end

export smatrix_real_space_tdvp

"""
    smatrix_expansion_real_space_tdvp(d, H, ψ0, w, χmax, t) -> Array{ComplexF64,4}

Prototype helper that computes the TDVP overlap tensor `A[j1, j2, j1′, j2′]` for two-particle states.

This is currently used as an intermediate object for momentum-space element evaluation.

# Warnings
- Experimental / under construction.
"""
function smatrix_expansion_real_space_tdvp(
    d::dType, # The local dimension of the system
    H::MPO, # The Hamiltonian of the system
    ψ0::ψ0Type, # The vacuum state
    w::wType, # The Wannier creation operator
    χmax::Int, # The maximum bond dimension
    t::tType # The time at which we want to compute the S-matrix element
    ) where {
        dType <: Integer,
        ψ0Type <: MPS,
        wType <: MPO,
        tType <: Real
        }

    # TODO assert that the local dimension of the MPS is the same as the one of the MPO
    # TODO assert that ψ0 is the groundstate of the Hamiltonian H

    Lψ0 = length(ψ0)
    l = length(w)
    L = Int64(Lψ0/2) - l + 1 # We define the length of the particle positions

    # We construct the single-particle basis wj
    println("Constructing the single-particle basis wj...")
    W = Array{MPO}(undef, L)
    Loffset = Int64(Lψ0/4)
    for j in 1:L
        W[j] = insert_local(j - 1 + Loffset, deepcopy(w), Lψ0 - Loffset - (j - 1) - l, d)
        substitute_siteinds!(W[j], ψ0)
    end

    # We construct the asimptotic in states
    println("Constructing the asimptotic in states...")
    ψin = Array{MPS}(undef, L, L)
    for j1 in 1:L
        for j2 in 1:L
            ψtoadd = deepcopy(ψ0)
            ψtoadd = apply(W[j1], ψtoadd; maxdim = χmax)
            ψtoadd = apply(W[j2], ψtoadd; maxdim = χmax)
            normalize!(ψtoadd)
            ψin[j1, j2] = ψtoadd
        end
    end

    # We construct the asimptotic out states as a deepcopy of the in states
    # We need also to substitute the site indices of the out states
    # println("Constructing the asimptotic out states...")
    ψout = deepcopy(ψin)
    # It seems that the following code is not needed...
    # for j1 in 1:L
    #     for j2 in 1:L
    #         substitute_siteinds!(ψout[j1, j2], ψ0)
    #     end
    # end

    # We compute the matrix A
    println("Computing the smatrix element expansion in real space...")
    A = zeros(ComplexF64, L, L, L, L) # N+1 because we include the N = 0 term
    for j1 in 1:L
        for j2 in 1:L
            if abs(j1 - j2) < l
                continue
            end
            println("TDVP Step: $j1, $j2")
            ψin[j1, j2] = tdvp(H, -im * t, ψin[j1, j2]; maxdim = χmax)
            for j1′ in 1:L
                for j2′ in 1:L
                    if abs(j1′ - j2′) < l
                        continue
                    end
                    A[j1, j2, j1′, j2′] = inner(ψout[j1′, j2′], ψin[j1, j2])
                end
            end
        end
    end

    return A
end

# Given A with 4 j indices and 1 n index, we can compute an S-matrix element in the momentum space
"""
    smatrix_element_momentum_space_tdvp(A, k1, k2, k1′, k2′) -> Complex

Compute a momentum-space scattering amplitude from a TDVP overlap tensor `A`.

# Arguments
- `A`: A 4-index overlap tensor (typically output of `smatrix_expansion_real_space_tdvp`).
- `k1, k2, k1′, k2′`: Momenta (radians).

# Returns
- A complex amplitude.

# Notes
- This uses an unnormalized direct summation; depending on conventions you may need to apply
  normalization factors externally.
"""
function smatrix_element_momentum_space_tdvp(A, k1, k2, k1′, k2′)
    result = 0
    L = size(A, 1)
    for j1 in 1:L
        for j2 in 1:L
            for j1′ in 1:L
                for j2′ in 1:L
                    result += exp(-im * k1 * j1) * exp(-im * k2 * j2) * exp(im * k1′* j1′) * exp(im * k2′ * j2′) * A[j1, j2, j1′, j2′]
                end
            end
        end
    end
    return result
end

export smatrix_element_momentum_space_tdvp

"""
    two_particle_wavefunction(ψ, ψ0, W, d; maxdim=100) -> Matrix{ComplexF64}

Compute a two-particle real-space wavefunction matrix from an MPS.

For each pair of insertion positions `(j, j′)`, this routine applies the creation operator `W` to the
vacuum `ψ0` at `j` and `j′` and computes the overlap with `ψ`.

# Arguments
- `ψ::MPS`: Target state (typically a two-particle state).
- `ψ0::MPS`: Vacuum/reference state.
- `W::MPO`: Creation operator MPO with support `length(W)`.
- `d::Int`: Local dimension.

# Keyword Arguments
- `maxdim::Int=100`: Maximum bond dimension used during `apply`.

# Returns
- A complex matrix `M[j, j′]` of overlaps.

# Warnings
- This routine is computationally expensive (`O(L^2)`) state preparations/overlaps.
"""
function two_particle_wavefunction(ψ::MPS, ψ0::MPS, W::MPO, d::Int; maxdim::Int = 100)
    normalize!(ψ)
    normalize!(ψ0)
    Lψ = length(ψ)
    L0 = length(ψ0)
    if siteinds(ψ) != siteinds(ψ0)
        substitute_siteinds!(ψ, ψ0)
    end
    @assert Lψ == L0
    Lw = length(W)
    function Wcr(j)
        ret = insert_local(j - 1, deepcopy(W), L0 - (j - 1) - Lw, d)
        substitute_siteinds!(ret, ψ0)
        return ret
    end
    L = Lψ - Lw + 1
    M = zeros(ComplexF64, L, L)
    for j in 1:L
        println("Step: $j")
        ψ0′ = deepcopy(ψ0)
        Wj = Wcr(j)
        ψ0′ = apply(Wj, ψ0′; maxdim = maxdim)
        normalize!(ψ0′)
        for j′ in 1:L
            if abs(j - j′) < Lw
                continue
            end
            ψ0′′ = deepcopy(ψ0′)
            ψ0′′ = apply(Wcr(j′), ψ0′′; maxdim = maxdim)
            normalize!(ψ0′′)
            M[j, j′] = inner(ψ0′′, ψ)
        end
    end
    return M
end
