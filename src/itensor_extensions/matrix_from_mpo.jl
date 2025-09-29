
"""
    matrix_from_mpo(mpo) -> Matrix

Converts an `MPO` to a `Matrix`, keeping the order of the indices defined in the `MPO`.

# Example
    julia> sx = [0 1; 1 0]
    julia> sz = [1 0; 0 -1]
    julia> mpo = mpo_from_matrix(kron(sx, sz, sx), 2)
    julia> matrix = matrix_from_mpo(mpo)
    8×8 Matrix{Float64}:
    0.0  0.0   0.0   0.0  0.0  1.0   0.0   0.0
    0.0  0.0   0.0   0.0  1.0  0.0   0.0   0.0
    0.0  0.0   0.0   0.0  0.0  0.0   0.0  -1.0
    0.0  0.0   0.0   0.0  0.0  0.0  -1.0   0.0
    0.0  1.0   0.0   0.0  0.0  0.0   0.0   0.0
    1.0  0.0   0.0   0.0  0.0  0.0   0.0   0.0
    0.0  0.0   0.0  -1.0  0.0  0.0   0.0   0.0
    0.0  0.0  -1.0   0.0  0.0  0.0   0.0   0.0
"""
function matrix_from_mpo(mpo::MPO)::Matrix
    L = length(mpo)
    sitesmain = siteinds_main(mpo)
    sitesprimed = siteinds_primed(mpo)
    tens = mpo[1]
    for j in 2:L
        tens = tens * mpo[j]
    end
    Cmain = combiner(sitesmain...)
    indmain = combinedind(Cmain)
    Cprimed = combiner(sitesprimed...)
    indprimed = combinedind(Cprimed)
    tens = tens * Cmain
    tens = tens * Cprimed
    return Matrix(tens, indprimed, indmain)
end

export matrix_from_mpo