# =============================================================================
# Spectrum.jl  --  finite-volume energy-momentum spectrum
# =============================================================================
#
# Turns the raw Bloch states (one per (k, level)) into a flat list of `Level`
# records carrying the integer momentum label
#
#     K = 2π Kint / L .
#
# The library returns only the reduced-zone half k ∈ [0, π] (koverpi ∈ [0,1]);
# since ε(-k)=ε(K) by reflection symmetry that is all we need for K = 0 work and
# for the one-particle band.  Negative-K sectors can be reconstructed by mirror.
# =============================================================================

# ----------------------------------------------------------------------------
# Output record (assignment-specified).  Kint means K = 2π Kint / L.
# ----------------------------------------------------------------------------
struct Level
    L::Int
    Kint::Int
    K::Float64
    index::Int        # 1 = lowest in this (L, Kint) sector, 2 = next, ...
    energy::Float64
end

"""
    compute_momentum_spectrum(H, T, L; nlevels, parity_op=nothing) -> Vector{Level}

Diagonalize and label every state by its lattice momentum.  Returns a flat,
energy-sorted-per-sector vector of `Level`.

The momentum integer is recovered from `koverpi = k/π = 2 Kint / L`, i.e.
`Kint = round(koverpi * L / 2)`.

If `parity_op` (a spatial-reflection operator) is supplied, then in the K=0
sector only reflection-EVEN states are kept.  This isolates the interacting
scattering channel for models where opposite parities are different channels
(Hubbard: even = singlet sees U, odd = triplet is free).  The vacuum and the
k=0 one-particle state are reflection-even, so the band is unaffected.
"""
function compute_momentum_spectrum(H, T, L::Int; nlevels::Int, parity_op = nothing)
    states = raw_dispersion(H, T, L; nlevels = nlevels)

    # group by integer momentum
    levels = Level[]
    bysector = Dict{Int, Vector{Float64}}()
    for s in states
        Kint = round(Int, Float64(s.koverpi) * L / 2)
        if parity_op !== nothing && Kint == 0
            ψ = s.data
            par = real(dot(ψ, parity_op * ψ)) / real(dot(ψ, ψ))
            par < 0 && continue        # drop reflection-odd (free) states at K=0
        end
        push!(get!(bysector, Kint, Float64[]), Float64(s.energy))
    end
    for (Kint, energies) in bysector
        sort!(energies)
        K = 2π * Kint / L
        for (i, E) in enumerate(energies)
            push!(levels, Level(L, Kint, K, i, E))
        end
    end
    sort!(levels, by = lv -> (lv.Kint, lv.index))
    return levels
end

# ----------------------------------------------------------------------------
# Convenience selectors
# ----------------------------------------------------------------------------
levels_in_sector(levels::Vector{Level}, Kint::Int) =
    sort(filter(lv -> lv.Kint == Kint, levels), by = lv -> lv.energy)

momentum_sectors(levels::Vector{Level}) = sort(unique(lv.Kint for lv in levels))

"""
    vacuum_energy(levels) -> (E0, Kint0)

Identify the vacuum as the global lowest level and return its energy and the
sector it lives in.  Warns if the vacuum is not in the K=0 sector or if it is
near-degenerate with another state (signalling a non-unique / symmetry-broken
vacuum, e.g. the ordered phase).
"""
function vacuum_energy(levels::Vector{Level})
    @assert !isempty(levels) "empty spectrum"
    i0 = argmin(lv.energy for lv in levels)
    E0 = levels[i0].energy
    Kint0 = levels[i0].Kint

    if Kint0 != 0
        @warn "Vacuum is not in the K=0 sector (found Kint=$Kint0). Unusual for a gapped phase."
    end

    # near-degeneracy check across the whole spectrum
    others = sort([lv.energy for lv in levels if lv.energy > E0 + 1e-9])
    if !isempty(others)
        gap = others[1] - E0
        if gap < 1e-2
            @warn "Vacuum near-degenerate: gap to next state = $(round(gap, sigdigits=3)). " *
                  "Likely a symmetry-broken / ordered phase (second ferromagnetic vacuum). " *
                  "One-particle identification will be unreliable here."
        end
    end
    return E0, Kint0
end
