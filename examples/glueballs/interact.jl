using GLMakie

# Define the function to plot
f(x, k) = sin.(k .* x)

# Create a figure and axis
fig = Figure()
ax = Axis(fig[1, 1])

# Define the x-range
x = LinRange(0, 2π, 1000)

# Initialize the plot with k = 1
k = Observable(1.0)
lineplot = lines!(ax, x, f(x, k[]))

# Create a slider to adjust k
slider = Slider(fig[2, 1], range = 0.1:0.1:10.0, startvalue = k[])

# Update the plot when the slider value changes
on(slider.value) do new_k
    k[] = new_k
    lineplot[1][] = f(x, k[])
end

# Display the figure
fig