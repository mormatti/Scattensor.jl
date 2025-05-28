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

export fourier_transform

function fourier_transform(M::Matrix)
    return fourier_transform(M, size(M, 1), size(M, 2))
end

export fourier_transform