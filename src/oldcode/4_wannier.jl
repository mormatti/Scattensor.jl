using SymPy

# Define the symbolic imaginary unit
imma = SymPy.I

# Define the cyclic sub-diagonal matrix T
T = SymPy.Matrix([[0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 0]])

# Compute T^dagger (transpose for real entries)
T_dagger = T.transpose()

# Compute iT - iT^dagger
M = imma * T - imma * T_dagger

println("Matrix iT - iT^dagger:")
println(M)

# Diagonalize the matrix
P, D = M.diagonalize()

println("Eigenvector matrix P:")
println(P)

println("Diagonal matrix D (Eigenvalues):")
println(D)