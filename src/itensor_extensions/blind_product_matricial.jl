"""
    blind_product_matricial(mpo1, mpo2) -> MPO
    blind_product_matricial(mpo, mps) -> MPS

Computes the inner product of two `MPS` and `MPO` objects
of same length without keeping track of the site indices instances, but just 
of their local dimensions.

# Example
    julia> mps = random_mps(siteinds(3, 5))
    julia> mpo = random_mpo(siteinds(3, 5))
    julia> blind_product_matricial(mpo, mps)
    MPS
    [1] ((dim=1|id=833|"Link,l=1"), (dim=3|id=224|"Site,n=1"))
    [2] ((dim=3|id=922|"Site,n=2"), (dim=1|id=773|"Link,l=2"), (dim=1|id=833|"Link,l=1"))
    [3] ((dim=3|id=48|"Site,n=3"), (dim=1|id=339|"Link,l=3"), (dim=1|id=773|"Link,l=2"))
    [4] ((dim=3|id=882|"Site,n=4"), (dim=1|id=792|"Link,l=4"), (dim=1|id=339|"Link,l=3"))
    [5] ((dim=3|id=669|"Site,n=5"), (dim=1|id=792|"Link,l=4"))
"""
function blind_product_matricial(mpo::MPO, mps::MPS)
    substitute_siteinds!(mps, mpo)
    return apply(mpo, mps)
end

function blind_product_matricial(mpo1::MPO, mpo2::MPS)
    substitute_siteinds!(mpo1, mpo2)
    return apply(mpo1, mpo2)
end

export blind_product_matricial