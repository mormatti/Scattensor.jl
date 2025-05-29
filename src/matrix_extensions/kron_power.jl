"""
    kron_power(A::AbstractMatrix, n::Integer) -> AbstractMatrix

`kron_power` implementation for a generic AbstractMatrix type.
Returns a matrix of the same type of `A`.
If `n` is zero, it returns the identity matrix 1x1 of size (1) of the same type of `A`.

# Examples
    julia> kron_power([1 2; 3 4], 2)
    4×4 Matrix{Int64}:
    1  2   2   4
    3  4   6   8
    3  6   4   8
    9 12  12  16

    julia> kron_power([1.0 2.0; 3.0 4.0], 0)
    1×1 Matrix{Int64}:
    1
    """
function kron_power(A::AbstractMatrix, n::Integer)
    if n == 0
        return operator_identity(typeof(A), 1)  # Returns the identity matrix of size A
    elseif n < 0
        error("The kronecker power with negative exponent is not defined.")
    else
        result = A
        for _ in 2:n
            result = kron(result, A)
        end
        return result
    end
end