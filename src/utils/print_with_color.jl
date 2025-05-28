"""A function to print colored text in the standard output."""
function print_with_color(text::String, color::RGB)
    r = round(Int, red(color) * 255)
    g = round(Int, green(color) * 255)
    b = round(Int, blue(color) * 255)
    print("\e[1m\e[38;2;$r;$g;$b;249m", text)
end

export print_with_color