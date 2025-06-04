# TODO: improve this function, considering a method which does not implement compression.

"""
    summation_local(mpo, L; convolution = identity, kwargs...)

Creates a summation of an `MPO` with uniform local dimensions, extending it to a total length `L`.
The `convolution` function is applied to the indices of the MPO, and it must return real or complex numbers.
In other words, it computes the sum `ΣⱼcⱼAⱼ`, where `Aⱼ` are the MPO tensors and `cⱼ` are the coefficients returned by the `convolution` function.
If the `convolution` function is not provided, it defaults to the identity function, which returns 1 for each index.
"""
function summation_local(mpo::MPO, L::Integer; convolution::Function = identity, kwargs...)
    if !is_uniform_localdim(mpo)
        error("The MPO must have uniform local dimensions.")
    end
    L0 = length(mpo)
    Lc = L - L0
    # We identify the convolution function
    if convolution == identity
        f(j::Integer) = 1
        convolution = f
    else
        if !hasmethod(sin, Tuple{Int})
            error("The convolution function must accept an integer argument.")
        end
    end
    # We assert that the convolution function returns only complex or real numbers
    coefficients = []
    for j in 1:Lc
        el = convolution(j)
        if !(el <: Union{Real, Complex})
            error("The convolution function must return only real numbers.")
        end
        push!(coefficients, el)
    end
    finalsum = coefficients[1] * insert_local(1 - 1, mpo, Lc)
    for j in 2:(Lc+1)
        tosum = coefficients[j] * insert_local(j - 1, mpo, Lc - j + 1)
        substitute_siteinds!(tosum, finalsum)
        finalsum += tosum
    end
    finalsum = truncate(finalsum, kwargs...)
end

export summation_local