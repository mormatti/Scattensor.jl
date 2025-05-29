"""
    kron_power(A, n)

Compute the n-th Kronecker power `A ⊗ A ⊗ A ⊗ ... ⊗ A` (`n` times) of an object `A` which supports a Kronecker product operation `kron`.
Types which does not implement `kron` directly must be dispatched to a specific implementation.
"""
function kron_power(A, n::Integer)
    if n <= 0
        error("Function kron_power does not implement zero or negative powers specifically for $(typeof(A)). Consider to implement kron_power(A::$(typeof(A)), n::Integer).")
    elseif n < 0
        error("Negative powers of matrices are not supported in function kron_power.")
    else
        result = A
        for _ in 2:n
            result = kron(result, A)
        end
        return result
    end
end