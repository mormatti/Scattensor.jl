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
    if checklength && length(mpo1) != length(mpo2)
        error("The two MPOs must have the same length.")
    end
    si1 = siteinds(mpo1)
    si2 = siteinds(mpo2)
    for j in eachindex(mpo1)
        old = si1[j]
        new = si2[j]
        old == new && continue
        @assert length(old) == length(new) "At site $j, different number of site indices"
        # Pair old -> new (bra->bra, ket->ket) by position
        pairs = map(=>, old, new)
        mpo1[j] = replaceinds(mpo1[j], pairs...)
    end
    return mpo1
end

function substitute_siteinds!(mps1::MPS, mps2::MPS; checklength::Bool = true)
    if checklength && length(mps1) != length(mps2)
        error("The two MPS must have the same length.")
    end
    si1 = siteinds(mps1)
    si2 = siteinds(mps2)
    for j in eachindex(mps1)
        old = si1[j]
        new = si2[j]
        old == new && continue
        mps1[j] = replaceinds(mps1[j], old => new)
    end
    return mps1
end

# substitute_siteinds(mps1::MPS, mps2::MPS; kwargs...) =
#     substitute_siteinds!(copy(mps1), mps2; kwargs...)

function substitute_siteinds!(mps::MPS, mpo::MPO; checklength::Bool = true, which::Symbol = :ket)
    if checklength && length(mps) != length(mpo)
        error("The MPS and MPO must have the same length.")
    end
    si_mps = siteinds(mps)
    si_mpo = siteinds(mpo)

    sel = which === :ket ? (j->si_mpo[j][2]) :
          which === :bra ? (j->si_mpo[j][1]) :
          error("`which` must be :ket or :bra")

    for j in eachindex(mps)
        old = si_mps[j]
        new = sel(j)
        old == new && continue
        mps[j] = replaceinds(mps[j], old => new)
    end
    return mps
end

#=
function substitute_siteinds!(mpo::MPO, mps::MPS; checklength::Bool = true)
    if checklength
        if length(mpo) != length(mps)
            error("The MPO and the MPS must have the same length.")
        end
    end
    simpo = siteinds(mpo)
    simps = siteinds(mps)
    for j in eachindex(mpo)
        replaceind!(mpo[j], simpo[j][1], prime(simps[j]))
        replaceind!(mpo[j], simpo[j][2], simps[j])
    end
end
=#

# substitute_siteinds deprecated!
# Now use the extension replace_siteinds

function ITensorMPS.replace_siteinds!(mpo::MPO, sites::AbstractVector{<:Index}; checklength::Bool = true)
    if checklength && length(mpo) != length(sites)
        error("MPO length must match length(sites).")
    end
    si = siteinds(mpo)
    for j in eachindex(mpo)
        inds = si[j]
        @assert length(inds) == 2 "MPO at site $j must have two site indices"
        old_ket, old_bra = inds[2], inds[1]
        new_ket = sites[j]
        new_bra = prime(sites[j])
        pairs = Pair{Index,Index}[]
        old_ket != new_ket && push!(pairs, old_ket => new_ket)
        old_bra != new_bra && push!(pairs, old_bra => new_bra)
        isempty(pairs) && continue
        mpo[j] = replaceinds(mpo[j], pairs...)
    end
    return mpo
end