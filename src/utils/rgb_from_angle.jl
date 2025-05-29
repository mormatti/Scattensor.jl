"""
    angle_to_rgb(angle::Float64) -> RGB

Converts an angle in radians `angle` to an RGB color using the HSL color model.
The hue is based on the angle, with full saturation and medium lightness.

# Example
    julia> angle_to_rgb(π/4)
    RGB{Float64}(0.5, 0.5, 1.0)
"""
function angle_to_rgb(angle::Float64)::RGB
    # Normalize angle to [0, 2π)
    hue = mod(angle, 2π)

    # Convert hue from radians to [0, 360) degrees for HSL
    hue_deg = rad2deg(hue)

    # Create HSL color and convert to RGB
    hsl_color = HSL(hue_deg, 1.0, 0.5)
    return RGB(hsl_color)
end

export angle_to_rgb