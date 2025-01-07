mutable struct EmptyLocalOperator
end
export ExactDiagLocalOperator

local_operator() = EmptyLocalOperator()

mutable struct PositionedLocalOperator
    matrix    :: Matrix{ComplexF64}
    position  :: Int
end

local_operator(matrix::Matrix{ComplexF64}, position::Int) = PositionedLocalOperator(matrix, position)
matrix(𝒪::EmptyLocalOperator) = 𝒪.matrix
position(𝒪::EmptyLocalOperator) = 𝒪.position
size(𝒪::EmptyLocalOperator) = size(𝒪.matrix)

mutable struct ExactDiagLocalOperator
    parent_system             :: ExactDiagSystem
    positioned_local_operator :: PositionedLocalOperator
end
export ExactDiagLocalOperator

function local_operator(𝒮::ExactDiagSystem, 𝐋::PositionedLocalOperator)
    return ExactDiagLocalOperator(𝒮, local_operator(matrix, position))
end

function local_operator(𝒮::ExactDiagSystem, 𝐌::Matrix{ComplexF64}, j::Int)
    d = 𝒮.local_dimension
    L = 𝒮.system_size
    j = j ↻ L
    message = "The local operator must have the same dimension of the local space."
    @assert (size(𝐌) == (d, d)) message
    return ExactDiagLocalOperator(𝒮, 𝐌, j)
end

LocalOperator::Type = Union{EmptyLocalOperator, PositionedLocalOperator, ExactDiagLocalOperator}

function set_position(𝒪::LocalOperator, position::Int)
    L = 𝒪.parent_system.system_size
    𝒪.position = position ↻ L
end

function set_matrix(𝒪::ExactDiagLocalOperator, matrix::Matrix{ComplexF64})
    d = 𝒪.parent_system.local_dimension
    @assert (size(matrix) == (d, d)) "The local operator must have the same dimension of the local space."
    𝒪.matrix = matrix
end