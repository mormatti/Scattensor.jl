"""
    m ↻ n -> Real

Computes the periodic modulus of `m` with respect to `n`, where the result is always in the range from `1` to `n`.
Compatible also with real numbers, but `n` must be non-zero.

# Example
    julia> 0↻3
    3
    julia> 3↻3
    3
    julia> 4.2↻3
    1.2
    julia> -0.1↻3
    2.9
    """
function ↻(n::Real, m::Real)::Real
    if n == 0
        throw(ArgumentError("Variable n must be non-zero in function n ↻ m."))
    end
    modu = mod(n, m)
    if modu == 0
        return m
    else
        return modu
    end
end

export ↻

"""
    ZZ(n, p) -> Real

A shortcut notation for the equivalence p mod(n) but centered in zero.
The function `ZZ(n, p)` computes the periodic modulus of `p` with respect to `n` centered in zero.
`p` and `n` can also be real numbers (not only integers), but `n` must be non-zero.
The result is always in the range from `-n/2` to `n/2`.

# Example
    julia> ZZ(3, 0)
    0
    julia> ZZ(3, 1)
    1
    julia> ZZ(3, 3)
    0
    julia> ZZ(3, 1.6)
    -1.4
    """
function ZZ(n::Real, p::Real)::Real
    if n == 0
        throw(ArgumentError("Variable n must be non-zero in function ZZ(n, p)."))
    end
    ret = (2 * mod((2 * p + n) / 2, n) - n) / 2
    if isinteger(ret)
        return Int(ret)
    else
        return Float64(ret)
    end
end

export ZZ