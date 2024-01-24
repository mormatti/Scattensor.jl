"""Abstract type for a general model."""
abstract type Model end

"""Mutable struct for the Ising Model."""
mutable struct IsingModel <: Model
    spin_interaction    ::Real
    transverse_field    ::Real
    longitudinal_field  ::Real
end
export IsingModel

"""Mutable struct for the Heisenberg Model."""
mutable struct HeisenbergModel <: Model
    x_spin_interaction  ::Real
    y_spin_interaction  ::Real
    z_spin_interaction  ::Real
    x_field             ::Real
    y_field             ::Real
    z_field             ::Real
end
export HeisenbergModel