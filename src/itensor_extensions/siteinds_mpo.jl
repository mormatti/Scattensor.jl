# TODO: document these functions.

# Returns the non-primed indices of the MPO as a vector of vectors.
function siteinds_main(mpo::MPO)
    v = siteinds(mpo)
    return [[v[i][j] for i in eachindex(v)] for j in 1:length(v[1])][2]
end

# Returns the primed indices of the MPO as a vector of indices.
function siteinds_primed(mpo::MPO)
    v = siteinds(mpo)
    return [[v[i][j] for i in eachindex(v)] for j in 1:length(v[1])][1]
end