"""
        plot_complex_colormap(Z::Matrix{<:Complex}; max_opacity = 1.0)

    Plots a complex matrix `Z` as a color map where the hue represents the angle of the complex numbers and the brightness represents their magnitude.
    The color map is generated using the `angle_to_rgb` function, which maps angles to RGB colors.

    # Arguments
    - `Z::Matrix{<:Complex}`: A matrix of complex numbers to be visualized.
    - `max_opacity::Float64=1.0`: The maximum opacity for the color map. Default is 1.0 (fully opaque).

    # Example
        julia> Z = rand(Complex{Float64}, 10, 10)
        julia> plot_complex_colormap(Z)
    
    # Function
"""
function plot_complex_colormap(Z::Matrix{<:Complex}; max_opacity = 1.0)
    mags = abs.(Z)
    max_mag = maximum(mags)
    mags_norm = max_mag == 0 ? mags : mags ./ max_mag
    angles = angle.(Z)
    color_matrix = [RGBA(angle_to_rgb(θ), α * max_opacity) for (θ, α) in zip(angles, mags_norm)]
    color_image = reshape(color_matrix, size(Z)...)
    Plots.heatmap(color_image, aspect_ratio=:equal, axis=nothing, xlims=(1/2, size(Z, 2) + 1/2), ylims=(1/2, size(Z, 1) + 1/2))
end

export plot_complex_colormap