"""A shortcut binary notation for the periodic modulus."""
function ↻(n::Integer, m::Integer)::Integer
    n > 0 ? (n-1)%m + 1 : m + n%m
end
export ↻

"""A shortcut binary notation for the periodic modulus centered in zero"""
function ZZ(n::Integer, p::Integer)::Integer
    s = ((n+1)/2)
    p = p + s
    return (p > 0 ? (p-1)%n + 1 : n + p%n) - s
end
export ZZ

