function substitute_siteinds!(mpo1::MPO, mpo2::MPO)
    @assert length(mpo1) == length(mpo2) "The MPOs must have the same length."
    si1 = siteinds(mpo1)
    si2 = siteinds(mpo2)
    for j in eachindex(mpo1)
        mpo1[j] = mpo1[j] * delta(si1[j][1], si2[j][1])
        mpo1[j] = mpo1[j] * delta(si1[j][2], si2[j][2])
        # TODO: nonsense to do that (it affects performance), one should do sometinhg like:
        # ITensors.replaceind(mpo1[j], si1[j][1], si2[j][1])
        # ITensors.replaceind(mpo1[j], si1[j][2], si2[j][2])
        # But dunno it does not work! Fix it.
    end
end

function substitute_siteinds!(mps1::MPS, mps2::MPS)
    @assert length(mps1) == length(mps2) "The two MPS must have the same length."
    si1 = siteinds(mps1)
    si2 = siteinds(mps2)
    for j in eachindex(mps1)
        mps1[j] = mps1[j] * delta(si1[j], si2[j])
    end
end

function substitute_siteinds!(mps::MPS, mpo::MPO)
    error("Don't use this method, substitute mpo indices with the mps one instead.")
end

function substitute_siteinds!(mpo::MPO, mps::MPS)
    @assert length(mpo) == length(mps) "The MPS and the MPO must have the same length."
    simpo = siteinds(mpo)
    simps = siteinds(mps)
    for j in eachindex(mpo)
        mpo[j] = mpo[j] * delta(simpo[j][1], prime(simps[j]))
        mpo[j] = mpo[j] * delta(simpo[j][2], simps[j])
    end
end

export substitute_siteinds!