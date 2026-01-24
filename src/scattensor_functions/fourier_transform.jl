"""
    fourier_transform(M::Matrix, N1::Int, N2::Int) -> Matrix{ComplexF64}
    fourier_transform(M::Matrix) -> Matrix{ComplexF64}

Compute a (naive) 2D discrete Fourier transform of a matrix.

This routine implements the transform explicitly with nested loops and is intended for small
matrices or prototyping.

# Arguments
- `M::Matrix`: Input matrix of size `(L1, L2)`.
- `N1::Int`: Number of momentum points along the first dimension.
- `N2::Int`: Number of momentum points along the second dimension.

# Returns
- A `N1 × N2` complex matrix containing the transformed values, normalized by `1/sqrt(L1*L2)`.

# Notes
- This is an `O(N1*N2*L1*L2)` implementation.
- TODO in source: consider replacing with a fast FFT-based implementation when appropriate.
"""
function fourier_transform(M::Matrix, N1::Int, N2::Int)
    L1 = size(M, 1)
    L2 = size(M, 2)
    result = zeros(ComplexF64, N1, N2)
    for k1 in 1:N1
        for k2 in 1:N2
            p1 = (k1 - 1) * 2π / N1
            p2 = (k2 - 1) * 2π / N2
            for j in 1:L1
                for j′ in 1:L2
                    result[k1, k2] += M[j, j′] * exp(-im * p1 * j) * exp(-im * p2 * j′)
                end
            end
        end
    end
    return result * (1 / sqrt(L1 * L2))
end

"""
    fourier_transform(M::Matrix) -> Matrix{ComplexF64}

Convenience overload of [`fourier_transform(M, N1, N2)`](@ref) using `N1 = size(M,1)` and
`N2 = size(M,2)`.
"""
function fourier_transform(M::Matrix)
    return fourier_transform(M, size(M, 1), size(M, 2))
end
