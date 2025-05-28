"""
        print_with_color(text::String, color::RGB)

    A function to print colored text in the standard output.

    # Arguments
    - `text::String`: The text to be printed.
    - `color::RGB`: The color to be used for printing the text. The color should be an instance of the `RGB` type from the `ColorTypes` package.
    
    # Function
    """
function print_with_color(text::String, color::RGB)
    r = round(Int, red(color) * 255)
    g = round(Int, green(color) * 255)
    b = round(Int, blue(color) * 255)
    print("\e[1m\e[38;2;$r;$g;$b;249m", text)
end

export print_with_color