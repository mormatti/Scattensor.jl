# DEFINITION

mutable struct BlochState{T}
    repr::T # The state of the Bloch state
    energy::Real # The energy of the Bloch state
    kfraction::Rational # The momentum of the Bloch state (in Rational representation)
    # parityphase::Real
end

# PROPERTIES

wavefunction(bs::BlochState) = bs.repr
export wavefunction

energy(bs::BlochState) = bs.energy
export energy

momentum(bs::BlochState) = bs.kfraction * π
export momentum

# parityphase(bs::BlochState) = bs.parityphase