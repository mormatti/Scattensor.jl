"""
    apply_translation(mps; moveright = true, kwargs...)

Apply the translation of a finite MPS with uniform local dimension.
The translation is performed swapping couple of physical indices consecutively.
The keyword arguments are applied to the `swapbondsites` function of ITensor, so one
can pass cutoff for instance.

# Example
    julia> sites = siteinds(3,4)
    julia> psi = random_mps(sites)
    julia> phi = apply_translation(psi)
"""
function apply_translation(mps::MPS; moveright::Bool = true, kwargs...)
    if !is_uniform_localdim(mps)
        error("The input MPS must have have uniform local dimensions.")
    end
    L = length(mps)
    mpscopy = copy(mps)
    if moveright
        for j in L-1:-1:1
            mpscopy = swapbondsites(mpscopy, j, kwargs...)
        end
    else
        for j in 1:L-1
            mpscopy = swapbondsites(mpscopy, j, kwargs...)
        end
    end
    return mpscopy
end

export apply_translation