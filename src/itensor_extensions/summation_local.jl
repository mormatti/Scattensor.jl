# TODO: improve this function, considering a method which does not implement compression (hard task).

"""
    summation_local(mpo, L; convolution = identity, kwargs...)

Creates a summation of an `MPO` with uniform local dimensions, extending it to a total length `L`.
The `convolution` function is applied to the indices of the MPO, and it must return real or complex numbers.
In other words, it computes the sum `ΣⱼcⱼAⱼ`, where `Aⱼ` are the MPO tensors and `cⱼ` are the coefficients returned by the `convolution` function.
If the `convolution` function is not provided, it defaults to the identity function, which returns 1 for each index.
"""
function summation_local(mpo::MPO, L::Integer; convolution::Function = identity, pbc = false, cutoff = default_cutoff, maxdim = default_maxdim)
    if !is_uniform_localdim(mpo)
        error("The MPO must have uniform local dimensions.")
    end
    d = get_uniform_localdim(mpo)
    L0 = length(mpo)
    Lc = L - L0
    # If the concolution is the identity function, we define it as a function returning 1
    if convolution == identity
        f(j::Integer) = 1
        convolution = f
    end
    # We construct the coefficients list
    coefficients = []
    for j in 1:(Lc+1)
        el = convolution(j)
        if !(typeof(el) <: Union{Real, Complex})
            error("The convolution function must return only real or complex numbers.")
        end
        push!(coefficients, el)
    end
    println()
    print("step: 1")
    # We compute the final MPO summation
    finalsum = coefficients[1] * insert_local(0, mpo, Lc)
    finalsites = siteinds_main(finalsum)
    for j in 2:(Lc+1)
        print(", $j ")
        tosum = coefficients[j] * insert_local(j - 1, mpo, Lc - (j - 1))
        replace_siteinds!(tosum, finalsites)
        finalsum = add(finalsum, tosum, cutoff = cutoff, maxdim = maxdim)
    end
    println(".")
    # We finally sum also pbc contribution if needed
    if pbc
        println("and PBC contributions...")
        translationop = operator_translation(MPO, d, L)
        replace_siteinds!(translationop, finalsites)
        finalterm = coefficients[Lc+1] * insert_local(Lc, mpo, 0)
        replace_siteinds!(finalterm, finalsites)
        for _ in (Lc+2):L
            finalterm = product(adjoint_mpo(translationop), finalterm, translationop)
            finalsum = add(finalsum, finalterm, cutoff = cutoff, maxdim = maxdim)
        end
    end
    return finalsum
end
