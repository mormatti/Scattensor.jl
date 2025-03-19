# # DEFINITION

# """ The object which represents a local operator acting on a chain of sites.
#     It can be either a matrix or a MPO (or someting else in the future maybe).
    
#     ## Fields
#     - `repr` is the representation of the local operator, which can be a matrix or a MPO;
#     - `sitedims` is a tuple of integers representing the sites where the local operator acts.
#        For instance, the local Hamiltonian of the Ising model has sites (2,2).
#     """
# mutable struct LocalOperator{T} <: QuantumOperator{T}
#     repr::T # The representation of the local operator
#     localdims:: Vector{Int64} # The local dimensions of the sites
#     name::String # The symbol representing the name of the local operator
# end

# # PROPERTIES

# function Base.length(lo::T) where {T <: LocalOperator}
#     return length(lo.localdims)
# end

# function hasuniformlocaldims(lo::T) where {T <: LocalOperator}
#     return all(x -> x == lo.localdims[1], lo.localdims)
# end

# function representation(lo::T) where {T <: LocalOperator}
#     return lo.repr
# end

# # CONSTRUCTORS

# function LocalOperator(matrix::T, localdims::Vector{U}, name::String) where {T <: AbstractMatrix, U <: Unsigned}
#     # We first assert that the matrix is squared and we get the dimension
#     N, N2 = size(matrix)
#     @assert N == N2 "The matrix must be squared."

#     # We assert that N = product of the integers in the localdims
#     @assert N == prod(localdims) "The matrix size must be the product of the local dimensions."

#     return LocalOperator{T}(matrix, localdims, name)
# end

# function LocalOperator(mpo::MPO, name::String)
#     # We extract the site indices of the MPO
#     sind = siteinds(mpo)
#     L = length(mpo)

#     # We assert that the upper and lower indices of the MPO have the same dimension
#     for i in 1:L
#         if dim(sind[i][1]) != dim(sind[i][2])
#             error("Dimension of the upper and lower indices of MPO must be the same.")
#         end
#     end

#     # Now we extract the local dimensions
#     localdims = Tuple(dim(s[1]) for s in sites)

#     return LocalOperator{MPO}(mpo, localdims, name)
# end