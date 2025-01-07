# Tensor Networks functions

    """
        Plots the dispersion relation of the system.
        Inputs:
        - `k` is the `Vector` of momenta of the eigenstates;
        - `ℰ` is the list of energies of the eigenstates;
        - `filename` is the name of the file (with extension) where the plot will be saved;
        - `aspectRatio` is the ratio between the height and the width of the plot (default 1:1).
        """
    function plot_dispersion_relation(
        k::Vector,
        ℰ::Vector;
        file_name::String = "plot.png",
        aspect_ratio = 1.0,
        points_color = RGBA{Float64}(1, 0, 0, 1),
        background_color = RGBA{Float64}(1, 1, 1, 1),
        marker_size = 2,
        dpi = 300
        )

        @debug "Plotting the dispersion relation..."

        # Plotting the dispersion relation
        Plots.plot(
            k, 
            ℰ, 
            seriestype = :scatter, 
            markersize = marker_size, 
            legend = false, 
            xlabel = "k", 
            ylabel = "ℰ",
            color = points_color,
            background_color = background_color,
            dpi = dpi
            )
        Plots.plot!(size=(170,170 * aspect_ratio), dpi = dpi) # Setting the size of the plot
        Plots.savefig(file_name) # Saving the plot in a file
    end