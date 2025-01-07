# In this snippet, we define the Hamiltonian of the Ising model with transverse and longitudinal 
# fields using ITensors.

using ITensors

modelName = "Ising"

# We define the local Hilbert space in ITensor
ITensors.space(::SiteType$modelName) = 2

# The identity operator
ITensors.op(::OpName"Id",::SiteType"Ising") = [1 0; 0 1]

# The Pauli matrices
ITensors.op(::OpName"σx",::SiteType"Ising") = [0 1; 1 0]
ITensors.op(::OpName"σy",::SiteType"Ising") = [0 -im; im 0]
ITensors.op(::OpName"σz",::SiteType"Ising") = [1 0; 0 -1]

# The raising and lowering operators
ITensors.op(::OpName"σ+",::SiteType"Ising") = [0 1; 0 0]
ITensors.op(::OpName"σ-",::SiteType"Ising") = [0 0; 1 0]

# The Hamiltonian


