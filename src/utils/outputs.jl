"""A function to print colored text in the standard output."""
function printcolored(r, g, b, text)
    print("\e[1m\e[38;2;$r;$g;$b;249m", text)
end

export printcolored

function mathematica_format(matrix::AbstractMatrix)
    formatted_rows = ["{" * join(row, ", ") * "}" for row in eachrow(matrix)]
    return "{" * join(formatted_rows, ", ") * "}"
end
export mathematica_format

using Colors

""" Maps a complex number `z` to a color using HSV space:
    - Hue corresponds to the argument (angle) of `z`
    - Value (brightness) and Saturation depend on the magnitude of `z` (can be adjusted)

    # Arguments
    - `z`: A complex number

    # Returns
    - A color in RGB format
"""
function complex_to_color(z::Complex)
    θ = angle(z)  # Argument (angle) of complex number
    r = abs(z)    # Magnitude of complex number

    hue = mod(θ / (2π), 1)         # Normalize angle to [0, 1]
    saturation = 1.0               # Full saturation
    value = 1.0 - exp(-r)          # Map magnitude to value (you can adjust this for effect)

    return RGB(HSV(hue, saturation, value))
end
export complex_to_color