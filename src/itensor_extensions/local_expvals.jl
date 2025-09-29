"""
    local_expvals(mps, matrix; hermitian = true)
    local_expvals(mps, mpo; hermitian = true) -> Vector
    local_expvals(mpsvector, mpo; hermitian = true) -> Matrx

Calculates the local expectation values of an `MPS` with respect to a p-local operator given by an `MPO`.
If the local operator is the same length as the `MPS`, it returns the inner product of the `MPS` with the product of the `MPO` and the `MPS`.
If the local operator is shorter than the `MPS`, it computes the local expectation values for each position in the `MPS` where the operator can be applied.
If the local operator is longer than the `MPS`, it raises an error.
If the input is a vector of `MPS`, it computes the local expectation values for each `MPS` in the vector and returns a matrix of results.
"""
function local_expvals(mps::MPS, mpo::MPO; hermitian::Bool = true)
    # We check that mps and mpo have compatible siteinds
    if siteinds(mps) != siteinds_main(mpo)
        error("MPS and MPO must share the same site indices.")
    end
    # We get the local dimension of the MPS and the MPO and we check that they are uniform
    dmps = get_uniform_localdim(mps)
    dA0 = get_uniform_localdim(mpo)
    # We check that the local dimensions are the same
    if dmps != dA0
        error("Different local dimension of MPS (d = $dmps) and MPO (d = $dA0).")
    end
    d = dmps
    # We check that the local operator is not greater than the MPS
    Lc = length(mps) - length(mpo)
    if Lc < 0
        error("The length of the local operator cannot be greater than the one of the MPS.")
    # If the local operator is the same length as the MPS, we can just return the inner product
    elseif Lc == 0
        val = inner(mps', mpo, mps)
        val = hermitian ? real(val) : val
        return val
    end
    vals = []
    for j in 1:(Lc+1)
        Aext = insert_local(j - 1, mpo, Lc - j)
        replace_siteinds!(Aext, siteinds(mps))
        val = inner(mps', Aext, mps)
        val = hermitian ? real(val) : val
        push!(vals, val)
    end
    return vals
end

function local_expvals(mps::MPS, matrix::Matrix; hermitian::Bool = true)
    mpo = mpo_from_matrix(matrix, get_uniform_localdim(mps))
    replace_siteinds!(mpo, siteinds(mps))
    return local_expvals(mps, mpo, hermitian = hermitian)
end

function local_expvals(mpsvector::Vector{MPS}, mpo::MPO; hermitian::Bool = true)::Matrix
    if !all(mps -> length(mps) == length(first(mpsvector)), mpsvector)
        error("All the MPS must have the same length in order to construct a matrix")
    end
    timelapsed = 0
    matrix = zeros(length(mpsvector), length(mpsvector[1]) - length(mpo))
    N = length(mpsvector)
    for n in 1:N
        print("Step $n / $N. Estimated remaining time: $(timelapsed * (N - n)).")
        timelapsed = @elapsed begin
        matrix[n,:] = local_expvals(mpsvector[n], mpo, hermitian = hermitian)
        end
        print("\r\u001b[2K")
    end
    return matrix
end

export local_expvals