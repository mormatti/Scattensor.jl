"""
    siteinds_main(mpo::MPO) -> Vector{Index}
    siteinds_primed(mpo::MPO) -> Vector{Index}

Return the ket (prime level 0) and bra (prime level 1) site indices of an `MPO`.

These helpers are useful for constructing compatible `MPS`/`MPO` objects and for performing
site-index substitutions.

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
    list = [[v[i][j] for i in eachindex(v)] for j in 1:length(v[1])][2]
    for ind in list
        if plev(ind) != 0
            error("MPO has wrong prime levels in site indices.")
        end
    end
    return list
end


function siteinds_primed(mpo::MPO)
    v = siteinds(mpo)
    list = [[v[i][j] for i in eachindex(v)] for j in 1:length(v[1])][1]
    for ind in list
        if plev(ind) != 1
            error("MPO has wrong prime levels in site indices.")
        end
    end
    return list
end

