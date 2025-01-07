mutable struct IsingModelParameters
    spin_interaction    :: Real
    transverse_field    :: Real
    longitudinal_field  :: Real
end
export IsingModelParameters

ising_model(; J::Real = 1, hx::Real = 0, hz::Real = 0) = IsingModelParameters(J, hx, hz)
interaction(ℳ::IsingModelParameters) = ℳ.spin_interaction
spin_interaction(ℳ::IsingModelParameters) = ℳ.spin_interaction
transverse_field(ℳ::IsingModelParameters) = ℳ.transverse_field
longitudinal_field(ℳ::IsingModelParameters) = ℳ.longitudinal_field

IsingModel::Type = Union{IsingModelParameters}


mutable struct HeisenbergModelParameters <: Model
    x_spin_interaction  :: Real
    y_spin_interaction  :: Real
    z_spin_interaction  :: Real
    x_field             :: Real
    y_field             :: Real
    z_field             :: Real
end
export HeisenbergModel

heisenberg_model(; Jx::Real = 1, Jy::Real = 1, Jz::Real = 1, hx::Real = 0, hy::Real = 0, hz::Real = 0) = HeisenbergModel(Jx, Jy, Jz, hx, hy, hz)
interaction(ℳ::HeisenbergModel) = (ℳ.x_spin_interaction, ℳ.y_spin_interaction, ℳ.z_spin_interaction)
spin_interaction(ℳ::HeisenbergModel) = interaction(ℳ)
x_spin_interaction(ℳ::HeisenbergModel) = ℳ.x_spin_interaction
y_spin_interaction(ℳ::HeisenbergModel) = ℳ.y_spin_interaction
z_spin_interaction(ℳ::HeisenbergModel) = ℳ.z_spin_interaction
field(ℳ::HeisenbergModel) = (ℳ.x_field, ℳ.y_field, ℳ.z_field)
x_field(ℳ::HeisenbergModel) = ℳ.x_field
y_field(ℳ::HeisenbergModel) = ℳ.y_field
z_field(ℳ::HeisenbergModel) = ℳ.z_field

HeisenbergModel::Type = Union{HeisenbergModelParameters}