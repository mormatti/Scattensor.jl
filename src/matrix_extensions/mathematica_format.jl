""" Given a matrix, returns the representation in the format of Wolfram mathematica.
This allows to output matrices from julia code and insert them as input in mathematica.

# Arguments
- `matrix`: the input matrix.
"""
function mathematica_format(matrix::AbstractMatrix)
    error("Mathematica_format method not implemented for this matrix type.")
end

function mathematica_format(matrix::Matrix{Float64})
    formatted_rows = ["{" * join(row, ", ") * "}" for row in eachrow(matrix)]
    return "{" * join(formatted_rows, ", ") * "}"
end

export mathematica_format