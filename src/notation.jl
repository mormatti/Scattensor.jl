using LinearAlgebra

"""A shortcut binary notation for the Kronecker product."""
function ⊗(A, B)
    return kron(A,B)
end

"""A shortcut binary notation for the periodic modulus."""
function ↻(n::Integer, m::Integer)::Integer
    n > 0 ? (n-1)%m + 1 : m + n%m
end