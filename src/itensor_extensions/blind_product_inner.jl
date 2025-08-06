"""
    blind_product_inner(mps1, mps2)

Computes the inner product of two `MPS` objects, `mps1` and `mps2`
of same length without keeping track of the site indices instances, but just 
of their local dimensions. This is performed by replacing the site indices of `mps2` 
with those of `mps1`

# Example
    julia> mps1 = random_mps(siteinds(3, 5))
    julia> mps2 = random_mps(siteinds(3, 5))
    julia> blind_product_inner(mps2, mps1)
    0.13671228015559006
"""
function blind_product_inner(mps1::MPS, mps2::MPS)
    replace_siteinds!(mps2, siteinds(mps1))
    return inner(mps1, mps2)
end

export blind_product_inner