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

function hex_to_rgb(hex::String; scale_255=false)
    rgb = parse(Colorant, hex)
    if scale_255
        R = round(Int, red(rgb) * 255)
        G = round(Int, green(rgb) * 255)
        B = round(Int, blue(rgb) * 255)
        return RGB(R, G, B)
    else
        return RGB(red(rgb), green(rgb), blue(rgb))
    end
end

export hex_to_rgb

"""
    angle_to_rgb(angle_rad::Float64) -> RGB

Converts an angle in radians to an RGB color using the HSL color model.
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

"""
    complex_colormap_plot(Z::AbstractMatrix; max_opacity=1.0)

Plots a matrix of complex numbers using hue for phase and opacity for magnitude.
"""
function complex_colormap_plot(Z::AbstractMatrix; max_opacity=1.0)
    mags = abs.(Z)
    max_mag = maximum(mags)
    mags_norm = max_mag == 0 ? mags : mags ./ max_mag

    angles = angle.(Z)
    color_matrix = [RGBA(angle_to_rgb(θ), α * max_opacity) for (θ, α) in zip(angles, mags_norm)]

    color_image = reshape(color_matrix, size(Z)...)

    Plots.heatmap(color_image, aspect_ratio=:equal, axis=nothing, border=:none)
end

export complex_colormap_plot