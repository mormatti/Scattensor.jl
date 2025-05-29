"""
    print_with_color(text::String, color::RGB)

A function to print a `String` `text` as colored text in the standard output.
The color is specified using an `RGB` type from the `ColorTypes` package.
"""
function print_with_color(text::String, color::RGB)
    r = round(Int, red(color) * 255)
    g = round(Int, green(color) * 255)
    b = round(Int, blue(color) * 255)
    print("\e[1m\e[38;2;$r;$g;$b;249m", text)
end

export print_with_color