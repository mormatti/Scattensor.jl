"""
    insert_local(L1, mpo, L2) -> MPO

Embed a local MPO into a larger chain by padding with identity MPOs.

This is a convenience wrapper around `kron` that returns:
- `mpo` if `L1 == 0` and `L2 == 0`
- `kron(mpo, I(L2))` if `L1 == 0`
- `kron(I(L1), mpo)` if `L2 == 0`
- `kron(I(L1), kron(mpo, I(L2)))` otherwise

where `I(L)` is the identity MPO of length `L` with the same local dimension as `mpo`.

# Example
    julia> sites = siteinds(3,4)
    julia> localmpo = random_mpo(sites)
    julia> result_mpo = insert_local(4, localmpo, 2)
"""
function insert_local(L1::Integer, mpo::MPO, L2::Integer)
    # We get the local dimension of the MPO, returning an error if it is not uniform
    d = get_uniform_localdim(mpo)
    # We define a shortcut for the identity operator
    idm(L) = operator_identity(MPO, d, L)
    # We insert the local MPO at position L1 and L2
    if L1 == 0 && L2 == 0
        return mpo
    elseif L1 == 0
        return kron(mpo, idm(L2))
    elseif L2 == 0
        return kron(idm(L1), mpo)
    else
        return kron(idm(L1), kron(mpo, idm(L2)))
    end
end 
