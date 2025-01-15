# DEFINITION

mutable struct BlochState{T}
    repr::T # The state of the Bloch state
    energy::Real # The energy of the Bloch state
    momentum::Rational # The momentum of the Bloch state (in Rational representation)
end

# PROPERTIES

energy(bs::BlochState) = bs.energy

momentum(bs::BlochState) = bs.momentum