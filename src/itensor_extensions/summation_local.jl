"""
    summation_local(mpo::MPO, L::Integer; convolution=identity, pbc=false, cutoff=default_cutoff, maxdim=default_maxdim) -> MPO

Sum translated copies of a local MPO along a 1D chain.

Given a local MPO `mpo` of length `L0`, this function builds an MPO of total length `L` that
represents a linear combination of all translated placements of `mpo` within the length-`L` chain.
Coefficients are provided by `convolution(j)` for the placement starting at position `j`.

Conceptually, it computes `sum_j c_j * A_j`, where `A_j` is `mpo` inserted at position `j`.

# Arguments
- `mpo::MPO`: Local MPO to translate and sum.
- `L::Integer`: Target total chain length.

# Keyword Arguments
- `convolution=identity`: Function `j -> c_j` returning a real/complex coefficient for each placement.
  If left as `identity`, it is treated as the constant function `1`.
- `pbc::Bool=false`: If `true`, also includes periodic boundary contributions by translating the last placement
  using the translation MPO.
- `cutoff=default_cutoff`, `maxdim=default_maxdim`: Compression parameters passed to `add`.

# Returns
- An `MPO` of length `L`.

# Notes
- This routine prints a small progress indicator while summing terms.
- TODO in source: consider alternative constructions that reduce compression overhead.
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
