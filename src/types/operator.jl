"""Insert the description of the struct."""

mutable struct EmptyOperator
end

operator() = EmptyOperator()

mutable struct OperatorProperties
    hilbert_space_dimension :: Integer
    hermiticity             :: Bool
    unitarity               :: Bool
end

operator_properties() = OperatorProperties(false, false, 0)
operator_properties(::EmptyOperator) = operator_properties()
function operator_properties(𝐌::Matrix)
    hermiticity = (𝐌 == 𝐌')
    unitarity = (𝐌 * 𝐌' == 𝐌' * 𝐌 == I)
    hilbert_space_dimension = size(𝐌)[1]
    return OperatorProperties(hermiticity, unitarity, hilbert_space_dimension)
end
hermiticity(𝒪::OperatorProperties) = 𝒪.hermiticity
unitarity(𝒪::OperatorProperties) = 𝒪.unitarity
hilbert_space_dimension(𝒪::OperatorProperties) = 𝒪.hilbert_space_dimension

mutable struct MatrixFormOperator
    matrix              :: Matrix{ComplexF64}
    operator_properties :: OperatorProperties
end

operator(𝐌::Matrix{ComplexF64}) = MatrixFormOperator(𝐌, operator_properties(𝐌))
matrix(𝒪::MatrixFormOperator) = 𝒪.matrix
operator_properties(𝒪::MatrixFormOperator) = 𝒪.operator_properties
hermiticity(𝒪::MatrixFormOperator) = hermiticity(operator_properties(𝒪))
unitarity(𝒪::MatrixFormOperator) = unitarity(operator_properties(𝒪))

mutable struct ExactDiagMatrixFormOperator
    matrix_form_operator     :: MatrixFormOperator
    parent_system            :: ExactDiagSystem
    translational_invariance :: Bool
end
export ExactDiagOperator

operator(𝐌::MatrixFormOperator, 𝒮::ExactDiagSystem) = ExactDiagMatrixFormOperator(𝐌, 𝒮, translational_invariance(𝐌, 𝒮))
matrix(𝒪::ExactDiagMatrixFormOperator) = matrix(𝒪.matrix_form_operator)
operator_properties(𝒪::ExactDiagMatrixFormOperator) = operator_properties(𝒪.matrix_form_operator)
hermiticity(𝒪::ExactDiagMatrixFormOperator) = hermiticity(𝒪.matrix_form_operator)
unitarity(𝒪::ExactDiagMatrixFormOperator) = unitarity(𝒪.matrix_form_operator)
translational_invariance(𝒪::ExactDiagMatrixFormOperator) = 𝒪.translational_invariance

function translational_invariance(𝐌::Matrix{ComplexF64}, 𝒮::ExactDiagSystem)::Bool
    𝐓 = matrix(translation_operator(𝒮))
    return 𝐌 * 𝐓 == 𝐓 * 𝐌
end

"""
Generates the translation operator T for a chain of L sites with local dimension d.
The system is assumed to be uniform, i.e. the local dimension is the same for all sites.
The system, in order to perform a translation, must be in periodic boundary conditions.

Inputs:
- `L` is the number of sites of the chain.
- `d` is the local dimension.

Outputs:
- The translation operator `T`.
"""

translation_operator(d::Integer,L::Integer)::Matrix{ComplexF64}
    N::Integer = d^L
    𝐓::Matrix{ComplexF64} = zeros(N,N)

    Lst = []
    c = 0
    for _ in 1:d
        lst = []
        for _ in 1:(N/d)
            c = c + 1
            push!(lst, c)
        end
        push!(Lst, lst)
    end

    for indL in eachindex(Lst)
        lst = Lst[indL]
        for ind in eachindex(lst)
            j = lst[ind]
            𝐓[j, ((d*(j-1)+1)%N) + indL - 1] = 1
        end
    end

    return 𝐓
end