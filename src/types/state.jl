"""Insert the description of the struct."""

# TYPE DEFINITION

mutable struct ExactDiagState
    vector                  :: Vector{ComplexF64}
    parent_system           :: Union{ExactDiagSystem, Missing}
    average_energy          :: Union{Float64, Missing}
    average_momentum        :: Union{Float64, Missing}
end
export ExactDiagState