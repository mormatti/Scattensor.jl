# DEFINITION

mutable struct BlochState{T}
    repr::T # The state of the Bloch state
    energy::Real # The energy of the Bloch state
    kfraction::Rational # The momentum of the Bloch state (in Rational representation)
end

# PROPERTIES

wavefunction(bs::BlochState) = bs.repr

energy(bs::BlochState) = bs.energy

momentum(bs::BlochState) = bs.kfraction * π