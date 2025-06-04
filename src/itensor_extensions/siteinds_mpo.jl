"""
    siteinds_main(mpo) -> Vector{Vector{Index}}
    siteinds_primed(mpo) -> Vector{Vector{Index}}

Returns the non-primed and primed indices of an `MPO` as vectors of vectors of `Index`.
This allows for easy access to the site indices of the MPO.

# Example
    julia> mpo = random_mpo(siteinds(3, 5))
    julia> siteinds_main(mpo)
    5-element Vector{ITensors.Index{Int64}}:
    (dim=3|id=25|"Site,n=1")
    (dim=3|id=255|"Site,n=2")
    (dim=3|id=208|"Site,n=3")
    (dim=3|id=864|"Site,n=4")
    (dim=3|id=515|"Site,n=5")
    julia> siteinds_primed(mpo)
    5-element Vector{ITensors.Index{Int64}}:
    (dim=3|id=25|"Site,n=1")'
    (dim=3|id=255|"Site,n=2")'
    (dim=3|id=208|"Site,n=3")'
    (dim=3|id=864|"Site,n=4")'
    (dim=3|id=515|"Site,n=5")'
"""
function siteinds_main(mpo::MPO)
    v = siteinds(mpo)
    return [[v[i][j] for i in eachindex(v)] for j in 1:length(v[1])][2]
end

function siteinds_primed(mpo::MPO)
    v = siteinds(mpo)
    return [[v[i][j] for i in eachindex(v)] for j in 1:length(v[1])][1]
end