"""
    operator_translation(MPO, d, L) -> MPO

Generate the translation operator in MPO representation, for a system of `L` sites with local dimension `d`.
The result bond dimension is `d^2` since the swapping operator is two-site.

# Example
    julia> operator_translation(MPO, 2, 3)
    MPO
    [1] ((dim=2|id=13|"Site,n=1")', (dim=2|id=13|"Site,n=1"), (dim=4|id=855|"Link,l=1"))
    [2] ((dim=2|id=403|"Site,n=2")', (dim=2|id=403|"Site,n=2"), (dim=4|id=855|"Link,l=1"), (dim=4|id=879|"Link,l=2"))
    [3] ((dim=2|id=617|"Site,n=3")', (dim=2|id=617|"Site,n=3"), (dim=4|id=879|"Link,l=2"))
"""
function operator_translation(::Type{MPO}, d::Integer, L::Integer)
    # First tensor (A)
    at = Index(d)
    ab = Index(d)
    at′ = prime(at)
    ab′ = prime(ab)
    A = delta(at, at′) * delta(ab, ab′)
    combA = combiner(ab′, at′)
    ar = combinedind(combA)
    A = A * combA

    # Bulk tensor (B)
    bt = Index(d)
    bb = Index(d)
    bm = Index(d)
    bt′ = prime(bt)
    bb′ = prime(bb)
    bm′ = prime(bm)
    B = delta(bt, bt′) * delta(bb, bb′) * delta(bm, bm′)
    combBl = combiner(bm, bb′)
    bl = combinedind(combBl)
    combBr = combiner(bm′, bt′)
    br = combinedind(combBr)
    B = B * combBl
    B = B * combBr

    # Last tensor (C)
    ct = Index(d)
    cb = Index(d)
    ct′ = prime(ct)
    cb′ = prime(cb)
    C = delta(ct, ct′) * delta(cb, cb′)
    combC = combiner(ct′, cb′)
    cl = combinedind(combC)
    C = C * combC

    # We create the list of tensors for the MPO and the list of sites
    vectorlist = Vector{ITensor}(undef, L)
    sites = siteinds(d, L)
    links = [Index(d * d, "Link,l=$j") for j in 1:L-1]

    # Populate the frst tensor
    W = copy(A)
    W = W * delta(at, prime(sites[1]))
    W = W * delta(ab, sites[1])
    W = W * delta(ar, links[1])
    vectorlist[1] = W

    # Populate the bulk tensors
    for j in 2:L-1
        W = copy(B)
        W = W * delta(bt, prime(sites[j]))
        W = W * delta(bb, sites[j])
        W = W * delta(bl, links[j-1])
        W = W * delta(br, links[j])
        vectorlist[j] = W
    end

    # Populate the last tensor
    W = copy(C)
    W = W * delta(ct, prime(sites[L]))
    W = W * delta(cb, sites[L])
    W = W * delta(cl, links[L-1])
    vectorlist[L] = W

    return adjoint_mpo(MPO(vectorlist))
end

export operator_translation