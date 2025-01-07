include("2_states.jl")

# We select the first 20 states in setstates which have minimum energy
# setstates = sort(collect(setstates), by = state -> energy(state))
# print the energies and momenta of the first 20 states
# for state in setstates[1:10]
# println("E = $(energy(state)), k = $(momentum(state))")
# end

