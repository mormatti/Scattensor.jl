"""
    complex_to_rgb(z::Complex) -> RGB

Maps a `Complex` number `z` to a color using HSV space where:
- Hue corresponds to the argument (angle) of `z`
- Value (brightness) and Saturation depend on the magnitude of `z` (can be adjusted)

# Example
    julia> complex_to_rgb(1 + 1im)
    RGB{Float64}(0.7071067811865475, 0.7071067811865475, 1.0)
"""
function complex_to_rgb(z::Complex)::RGB
    θ = angle(z)  # Argument (angle) of complex number
    r = abs(z)    # Magnitude of complex number
    hue = mod(θ / (2π), 1)         # Normalize angle to [0, 1]
    saturation = 1.0               # Full saturation
    value = 1.0 - exp(-r)          # Map magnitude to value (you can adjust this for effect)
    return RGB(HSV(hue, saturation, value))
end

export complex_to_rgb