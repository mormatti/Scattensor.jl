"""
    substitute_siteinds!(mpo1, mpo2; checklength = true) -> MPO
    substitute_siteinds!(mps1, mps2; checklength = true) -> MPS
    substitute_siteinds!(mps, mpo; checklength = true) -> MPS
    substitute_siteinds!(mpo, mps; checklength = true) -> MPO

Substitutes the site indices of one MPO or MPS with those of another.
This function is useful when you want to ensure that two MPOs or MPSs have compatible site indices, 
especially when performing operations like contraction or applying operators.

# Examples
    julia> mps1 = random_mps(siteinds(3, 5))
    julia> mps2 = random_mps(siteinds(3, 5))
    julia> siteinds(mps1) == siteinds(mps2)
    false
    julia> substitute_siteinds!(mps1, mps2)
    julia> siteinds(mps1) == siteinds(mps2)
    true
"""
function substitute_siteinds!(mpo1::MPO, mpo2::MPO; checklength::Bool = true)
    if checklength
        if length(mpo1) != length(mpo2)
            error("The two MPOs must have the same length.")
        end
    end
    si1 = siteinds(mpo1)
    si2 = siteinds(mpo2)
    if si1 != si2
        for j in eachindex(mpo1)
            replaceind!(mpo1[j], si1[j][1], si2[j][1])
            replaceind!(mpo1[j], si1[j][2], si2[j][2])
        end
    end
end

function substitute_siteinds!(mps1::MPS, mps2::MPS; checklength::Bool = true)
    if checklength
        if length(mps1) != length(mps2)
            error("The two MPS must have the same length.")
        end
    end
    si1 = siteinds(mps1)
    si2 = siteinds(mps2)
    if si1 != si2
        for j in eachindex(mps1)
            replaceind!(mps1[j], si1[j], si2[j])
        end
    end
end

function substitute_siteinds!(mps::MPS, mpo::MPO; checklength::Bool = true)
    if checklength
        if length(mps) != length(mpo)
            error("The MPS and the MPO must have the same length.")
        end
    end
    simps = siteinds(mps)
    simpo = siteinds_main(mpo)
    if simps != simpo
        for j in eachindex(mps)
            replaceind!(mps[j], simps[j], simpo[j])
        end
    end
end

function substitute_siteinds!(mpo::MPO, mps::MPS; checklength::Bool = true)
    @warn "mpo -> mps indices substitution is not recommended, since it could be ambiguous in its definition. Prefer using mps -> mpo instead."
    if checklength
        if length(mpo) != length(mps)
            error("The MPO and the MPS must have the same length.")
        end
    end
    simpo = siteinds(mpo)
    simps = siteinds(mps)
    for j in eachindex(mpo)
        replaceind!(mpo[j], simpo[j][1], simps[j])
        replaceind!(mpo[j], simpo[j][2], prime(simps[j]))
    end
end

export substitute_siteinds!