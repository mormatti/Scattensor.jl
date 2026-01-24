"""
    adjoint_mpo(A::MPO) -> MPO

Return the adjoint (dagger) of an MPO, with site-index priming normalized.

This helper computes `dag(A)` and then swaps prime levels so that the returned MPO uses
prime level 1 for bra indices and prime level 0 for ket indices (consistent with the
conventions used in the rest of this package).

# Arguments
- `A::MPO`: Input MPO.

# Returns
- An `MPO` representing the adjoint of `A`.
"""
function adjoint_mpo(A::MPO)::MPO
    Ad = dag(A)
    return swapprime(Ad, 0 => 1; tags = "Site")
end
