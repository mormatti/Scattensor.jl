"""
    mutable struct BlochState{T}

A structure representing a Bloch state in a quantum system.
The field `data::T` gives the representation of the Bloch state, in whatever format (e.g. `T` can be a `Vector` or `MPS`).
The `energy::Real` is the eigenvalue associated with the Bloch state, while `koverpi::Rational` is the crystal momentum in units of π.
The `koverpi` is computed so that it falls within the range [-1, 1).
Expressing the momentum as units of π is useful for periodic systems, where the momentum is typically defined as a Rational type.
"""
mutable struct BlochState{T}
    data::T # The state of the Bloch state, in whatever representation.
    energy::Real # The energy of the Bloch state
    koverpi::Real # The momentum of the Bloch state (/π so it can be represented as Rational)
end

export BlochState

"""
    wavefunction(bs::BlochState) -> Any

Returns the data representation of the Bloch state, such as its wavefunction or coefficients.
"""
wavefunction(bs::BlochState) = bs.data

export wavefunction

"""
    energy(bs::BlochState) -> Real

Returns the corresponding energy of the Bloch state.
"""
energy(bs::BlochState)::Real = bs.energy

export energy

"""
    momentum(bs::BlochState) -> Real

Returns the corresponding momentum `k` ∈ [-π,π) of the Bloch state.
"""
momentum(bs::BlochState)::Real = π * bs.koverpi

export momentum