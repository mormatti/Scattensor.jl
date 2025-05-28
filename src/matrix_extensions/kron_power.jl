function kron_power(A, n)
    result = A
    for _ in 2:n
        result = kron(result, A)
    end
    return result
end

export kron_power