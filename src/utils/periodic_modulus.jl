"""
        ↻(n::Integer, m::Integer) -> Integer

    A shortcut binary notation for the periodic modulus.
    
    # Example
        julia> 0↻3
        3
        julia> 1↻3
        1
        julia> 2↻3
        2
        julia> 3↻3
        3
        julia> 4↻3
        1
    
    # Function
    """
function ↻(n::Integer, m::Integer)::Integer
    n > 0 ? (n-1)%m + 1 : m + n%m
end
export ↻

"""
        ZZ(n::Integer, p::Integer) -> Integer

    A shortcut notation for the periodic modulus centered in zero.
    
    # Example
        julia> ZZ(3,0)
        0
        julia> ZZ(3,1)
        1
        julia> ZZ(3,2)
        -1
        julia> ZZ(3,3)
        0
    
    # Function
    """
function ZZ(n::Integer, p::Integer)::Integer
    s = ((n+1)/2)
    p = p + s
    return (p > 0 ? (p-1)%n + 1 : n + p%n) - s
end
export ZZ