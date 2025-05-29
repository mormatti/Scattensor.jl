"""
    mathematica_format(matrix::Matrix{<:Real}) -> String

Given a `Matrix`, returns the representation in the format of Wolfram mathematica.
This allows to output matrices from julia code and insert them as input in mathematica.

# Example
    julia> mathematica_format([1 2; 3 4])
    "{{1, 2}, {3, 4}}"
"""
function mathematica_format(matrix::Matrix{<:Real})
    formatted_rows = ["{" * join(row, ", ") * "}" for row in eachrow(matrix)]
    return "{" * join(formatted_rows, ", ") * "}"
end

export mathematica_format