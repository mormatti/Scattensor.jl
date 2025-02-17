"""A shortcut binary notation for the periodic modulus."""
function ↻(n::Integer, m::Integer)::Integer
    n > 0 ? (n-1)%m + 1 : m + n%m
end

"""A shortcut binary notation for the periodic modulus centered in zero"""
function ZZ(n::Integer, p::Integer)::Integer
    s = ((n+1)/2)
    p = p + s
    return (p > 0 ? (p-1)%n + 1 : n + p%n) - s
end

"""A function to print colored text in the standard output."""
function printcolored(r, g, b, text)
    print("\e[1m\e[38;2;$r;$g;$b;249m", text)
end