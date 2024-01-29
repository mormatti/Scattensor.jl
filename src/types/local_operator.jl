"""Insert the description of the struct."""

# TYPE DEFINITION

mutable struct LocalOperator
    parent_system::ExactDiagSystem
    matrix::Matrix{ComplexF64}
    position::Int
end
export LocalOperator


# CONSTRUCTORS

function LocalOperator(𝒮::ExactDiagSystem, matrix::Matrix{ComplexF64}, position::Int)
    d = 𝒮.local_dimension
    L = 𝒮.system_size
    position = position ↻ L
    @assert (size(matrix) == (d, d)) "The local operator must have the same dimension of the local space."
    return LocalOperator(𝒮, matrix, position)
end


# STRUCT FUNCTIONS

function set_position(𝒪::LocalOperator, position::Int)
    L = 𝒪.parent_system.system_size
    𝒪.position = position ↻ L
end

function set_matrix(𝒪::LocalOperator, matrix::Matrix{ComplexF64})
    d = 𝒪.parent_system.local_dimension
    @assert (size(matrix) == (d, d)) "The local operator must have the same dimension of the local space."
    𝒪.matrix = matrix
end