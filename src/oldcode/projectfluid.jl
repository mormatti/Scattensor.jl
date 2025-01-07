using ITensors
N = 10
s = siteinds(2,N)
chi = 4
psi = random_mps(s;linkdims=chi)


# The following function takes as an imput an integer and returns it in the 
# binary notation. The binary notation is a list of 0 and 1.
function binary(n::Int)
    if n == 0
        return [0]
    end
    b = []
    while n > 0
        pushfirst!(b, n % 2)
        n = div(n, 2)
    end
    return b
end

# The same previous function but it accepts a parameter N which is the number of
# bits we want to represent the integer.
function binary(n::Int, N::Int)
    b = binary(n)
    while length(b) < N
        pushfirst!(b, 1)
    end
    return b
end

# This function, given a list (the basis element) 
function element(b::Array{Int,1})
    el = [1,2,1,1,2,1,2,2,2,1]
    V = ITensor(1.)
    for j=1:N
      global V *= (psi[j]*state(s[j],el[j]))
    end
    v = scalar(V) 
end

# v is the element we wanted to obtain:
@show v