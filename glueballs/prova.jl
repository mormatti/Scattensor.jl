function C2(m...)
    N = length(m)
    s = 0
    s += 1//2 * sum([m[k]^2 for k in 1:N])
    s -= 1//(2*N) * (sum([m[k] for k in 1:N]))^2
    s += sum([m[k] * ((N+1)//2 - k) for k in 1:N])
end

N = 3

function C2sp(N,n)
    n = mod(n, N)
    return ((N+1)/(2*N))*(N-n)*(n)
end

E2(N,n,m) = C2sp(N,n) + C2sp(N,m) + C2sp(N,n-m)

TwoNE2mat = [E2(N,n,m) for n in 0:(N-1), m in 0:(N-1)] * N / (N + 1)

function bmat(gammamat)
    N = size(gammamat)[1]
    b = zeros(N, N)
    for n in 0:N-1
        in = n + 1
        in2 = mod(-n, N) + 1
        for m in 0:N-1
            im = m + 1
            im2 = mod(m - n, N) + 1
            b[in, im] = gammamat[in, im] * gammamat[in2, im2]
        end
    end
    return b
end
gamma2 = [-1 -1; +√(1/2) -√(1/2)]
gamma3 = [1 1 1; -√√(1/3) +√(1/3) -√√(1/3); +√(1/3) +√√(1/3) +√√(1/3)]
# gamma4 = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]
gamma5 = [  1 1 1 1 1; 
            1 1 1 1 1; 
            1 1 1 1 1; 
            1 1 1 1 1; 
            1 1 1 1 1]

bmat(gamma3)
