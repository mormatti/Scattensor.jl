""" Apply the translation of a finite MPS with uniform local dimension.
    The translation is performed swapping couple of physical indices consecutively.
    """
function apply_translation(ψ::MPS; dir = "right", cutoff = 1e-15)
    L = length(ψ)
    ϕ = copy(ψ)
    if dir == "left"
        for j in 1:L-1
            ϕ = swapbondsites(ϕ, j, cutoff = cutoff)
        end
    elseif dir == "right"
        for j in L-1:-1:1
            ϕ = swapbondsites(ϕ, j, cutoff = cutoff)
        end
    else
        error("Invalid direction")
    end
    return ϕ
end

export apply_translation