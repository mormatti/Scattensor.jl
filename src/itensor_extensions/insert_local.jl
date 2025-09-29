"""
    insert_local(L1, mpo, L2) -> MPO

Inserts a local MPO at positions L1 and L2 in a larger MPO, returning the resulting MPO.
The local MPO is inserted in the form of a direct product with identity operators at the other positions.

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

export insert_local