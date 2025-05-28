""" Converts an angle in radians to an RGB color using the HSL color model.
The hue is based on the angle, with full saturation and medium lightness.
"""
function angle_to_rgb(angle_rad::Float64)::RGB
    # Normalize angle to [0, 2π)
    hue = mod(angle_rad, 2π)

    # Convert hue from radians to [0, 360) degrees for HSL
    hue_deg = rad2deg(hue)

    # Create HSL color and convert to RGB
    hsl_color = HSL(hue_deg, 1.0, 0.5)
    return RGB(hsl_color)
end

export angle_to_rgb