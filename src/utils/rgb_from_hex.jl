"""
        hex_to_rgb(hex::String; scale_255 = false) -> RGB

    Convert a hexadecimal color string (e.g., `"#FF5733"`) to an `RGB` color object.

    # Arguments
    - `hex::String`: A string representing the hexadecimal color. It should be in the format `"#RRGGBB"` or `"RRGGBB"`.
    - `scale_255::Bool=false`: If `true`, the resulting RGB values will be scaled to the 0–255 integer range. If `false` (default), values will be in the 0.0–1.0 floating-point range.

    # Returns
    - An `RGB` object with red, green, and blue components either in the 0.0–1.0 or 0–255 range, depending on the `scale_255` keyword argument.

    # Examples
        julia> hex_to_rgb("#FF5733")
        RGB{Float64}(1.0, 0.3411764705882353, 0.2)

        julia> hex_to_rgb("#FF5733", scale_255=true)
        RGB{Int64}(255, 87, 51)
    """
function hex_to_rgb(hex::String; scale_255 = false)::RGB
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