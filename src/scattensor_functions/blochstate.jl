"""
    BlochState{T}

Container for a (single-particle) Bloch eigenstate.

`BlochState` stores a representation of the state (e.g. a coefficient vector or an `MPS`)
along with its energy and crystal momentum.

# Fields
- `data::T`: State representation (for example `Vector{ComplexF64}` or `ITensorMPS.MPS`).
- `energy::Real`: Energy eigenvalue associated with the state.
- `koverpi::Real`: Crystal momentum in units of π, i.e. `k/π`. In many workflows this is
  chosen in the reduced zone `[-1, 1)`, but this is a *convention* and is not enforced.

# Notes
- If you want exact momentum comparisons (e.g. grouping states by momentum), prefer storing
  `koverpi` as a `Rational` (which is a subtype of `Real`) rather than a floating-point number.
"""
mutable struct BlochState{T}
    data::T # The state of the Bloch state, in whatever representation.
    energy::Real # The energy of the Bloch state
    koverpi::Real # The momentum of the Bloch state (/π so it can be represented as Rational)
end

export BlochState

"""
    wavefunction(bs::BlochState) -> Any

Return the underlying representation stored in `bs.data` (e.g. a coefficient vector or an `MPS`).
"""
wavefunction(bs::BlochState) = bs.data

export wavefunction

"""
    energy(bs::BlochState) -> Real

Return the energy eigenvalue stored in `bs.energy`.
"""
energy(bs::BlochState)::Real = bs.energy

export energy

"""
    momentum(bs::BlochState) -> Real

Return the momentum `k` in radians, computed as `π * bs.koverpi`.

Note: no wrapping to a Brillouin zone is performed; the returned value follows whatever
convention was used to set `koverpi`.
"""
momentum(bs::BlochState)::Real = π * bs.koverpi

export momentum