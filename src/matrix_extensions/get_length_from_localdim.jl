# TODO: add documentation to these functions
function get_length_from_localdim(N::Integer, d::Integer)::Integer
    L = log(N) / log(d)
    # Check that L0 is an Integer
    if !(L ≈ round(L))
        error("Size of the local operator matrix incompatible with the local dimension d.")
    end
    # Round L to the nearest integer
    L = Int(round(L))
    return L
end

function get_length_from_localdim(vector::Vector, d::Integer)::Integer
    # Check if vector is a valid vector
    N = length(vector)
    if N == 0
        error("The input vector must not be empty.")
    end
    return get_length_from_localdim(N, d)
end

function get_length_from_localdim(matrix::Matrix, d::Integer)::Integer
    # Check if matrix is a square matrix
    N = size(matrix)[1]
    N2 = size(matrix)[2]
    if N != N2
        error("The input matrix must be square, i.e. have the same number of rows and columns.")
    end
    if N == 0
        error("The input matrix must not be empty.")
    end
    return get_length_from_localdim(N, d)
end