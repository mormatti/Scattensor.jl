"""
    hex_to_rgb(hex::String; scale_255 = false) -> RGB

Convert a hexadecimal color `hex` (e.g., `"#FF5733"`) to an `RGB` color object.
The format can be either `"#RRGGBB"` or `"RRGGBB"`.
By default, the RGB values are returned in the range of 0.0 to 1.0.
If `scale_255` is set to `true`, the RGB values will be scaled to the 0–255 integer range.

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