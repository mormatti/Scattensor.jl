"""
    plot_complex_colormap(Z::Matrix{<:Complex}; max_opacity = 1.0)
This another version of the function that takes a matrix of complex numbers and generates a color map based on the argument (angle) and magnitude of the complex numbers.
"""
function plot_complex_colormap(Z::Matrix{<:Complex})
    mags = abs.(Z)
    max_mag = maximum(mags)
    mags_norm = max_mag == 0 ? mags : mags ./ max_mag
    angles = angle.(Z)
    color_matrix = [RGBA(angle_to_rgb(θ), α * max_opacity) for (θ, α) in zip(angles, mags_norm)]
    color_image = reshape(color_matrix, size(Z)...)
    Plots.heatmap(color_image, aspect_ratio=:equal, axis=nothing, border=:none)
end

export plot_complex_colormap