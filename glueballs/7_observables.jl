# using Scattensor
# using ITensors, ITensorMPS
# using LinearAlgebra
# using SparseArrays
# using Optim
# using Plots
# using PlotlyJS
# using KrylovKit
# using Logging
# using Revise
# using JLD2
# plotlyjs()

# print("Loading data...")
# data = load("simulationdata.jld2")["data"]
# println("Done.")

# states = data[:states]

# print("Computing projectors...")
# projectors = Dict()
# projectors["121"] = Glueballs.projectr([1,2,1])
# projectors["131"] = Glueballs.projectr([1,3,1])
# projectors["1221"] = Glueballs.projectr([1,2,2,1])
# projectors["1231"] = Glueballs.projectr([1,2,3,1])
# projectors["1321"] = Glueballs.projectr([1,3,2,1])
# projectors["1331"] = Glueballs.projectr([1,3,3,1])
# print("Done.")

# observabledata = Dict()
# for numb in ("1221","1231","1321","1331")
#     println("Computing projector $numb")
#     observabledata[numb] = local_expvals(states, projectors[numb], 3)
# end
# print("Done.")

