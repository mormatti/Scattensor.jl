"""
    kron_power(A, n::Integer)

Compute the Kronecker power `A ⊗ A ⊗ ... ⊗ A` (n times).

This is the generic fallback that repeatedly applies `kron`. Concrete backends should provide
specialized methods when needed (e.g. to support `n == 0`, to handle identity objects, or to avoid
materializing large intermediates).

# Arguments
- `A`: Any object for which `kron(A, A)` is defined.
- `n::Integer`: Number of Kronecker factors.

# Returns
- The Kronecker power of `A` with itself `n` times.

# Notes
- The default method errors for `n <= 0`. Matrix and tensor backends may override this behavior.
"""
function kron_power(A, n::Integer)
    if n <= 0
        error("Function kron_power does not implement zero or negative powers specifically for $(typeof(A)). Consider to implement specific function kron_power(A::$(typeof(A)), n::Integer) with also n ≤ 0.")
    else
        result = A
        for _ in 2:n
            result = kron(result, A)
        end
        return result
    end
end
