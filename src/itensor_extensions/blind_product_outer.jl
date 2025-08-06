"""
    blind_product_outer(mps1, mps2) -> MPO

Computes the outer product |ψ1⟩⟨ψ2| of two `MPS` objects, `mps1` and `mps2`.
The final mpo will have the siteinds of `mps1`.

# Example
    julia> mps1 = random_mps(siteinds(3, 5))
    julia> mps2 = random_mps(siteinds(3, 5))
    julia> blind_product_outer(mps1, mps2)
    MPO
    [1] ((dim=3|id=662|"Site,n=1")', (dim=3|id=662|"Site,n=1"), (dim=1|id=71|"Link,l=1"))
    [2] ((dim=3|id=350|"Site,n=2")', (dim=3|id=350|"Site,n=2"), (dim=1|id=64|"Link,l=2"), (dim=1|id=71|"Link,l=1"))
    [3] ((dim=3|id=547|"Site,n=3")', (dim=3|id=547|"Site,n=3"), (dim=1|id=529|"Link,l=3"), (dim=1|id=64|"Link,l=2"))
    [4] ((dim=3|id=885|"Site,n=4")', (dim=3|id=885|"Site,n=4"), (dim=1|id=898|"Link,l=4"), (dim=1|id=529|"Link,l=3"))
    [5] ((dim=3|id=503|"Site,n=5")', (dim=3|id=503|"Site,n=5"), (dim=1|id=898|"Link,l=4"))
"""
function blind_product_outer(mps1::MPS, mps2::MPS)
    substitute_siteinds!(mps2, mps1)
    return outer(mps1', mps2)
end

export blind_product_outer