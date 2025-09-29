using SparseArrays
using Scattensor
using Revise
using Plots
using LinearAlgebra
using ITensorMPS, ITensors

# Global parameters
default_cutoff = 1e-10
default_maxdim = 100

# Use this
⊗(A, B) = kron(A, B)

# Ising model as mapping in ladder U(1) LGT
if false
    global d = 2
    global L0 = 3
    global Li = 7
    global λ = 0.5
    global id = SparseMatrixCSC([1 0; 0 1])
    global σx = SparseMatrixCSC([0 1; 1 0])
    global σz = SparseMatrixCSC([1 0; 0 -1])
    global hunit = 4 / π^2
    global ME = hunit * spdiagm(0 => [2, 2, 2, 0])
    global H0 = SparseMatrixCSC(λ * HE - (1 - λ) * HB)
    global Hi = summation_local(H0, d, Li; pbc = true)
    global Ti = operator_translation(SparseMatrixCSC, d, Li)
    global Pi = operator_reflection(SparseMatrixCSC, d, Li)
end

# Ising as deformation around critical point:
# H = (- ∑ σzσz - ∑ σx) - hx ∑ σx - hz ∑ σz
# Hloc = HZZ + (1+hx) HX + hz HZ
if false
    global d = 2
    global L0 = 3
    global Li = 17
    global hx = -0.5
    global hz = 0.5
    global id = SparseMatrixCSC([1 0; 0 1])
    global σx = SparseMatrixCSC([0 1; 1 0])
    global σz = SparseMatrixCSC([1 0; 0 -1])
    HZZ = -σz ⊗ σz
    HX = -0.5 * σx ⊗ id - 0.5 * id ⊗ σx
    HZ = -0.5 * σz ⊗ id - 0.5 * id ⊗ σz
    global Htwo = (HZZ + HX) + hx * HX + hz * HZ
    global H0 = 0.5 * (Htwo ⊗ id + id ⊗ Htwo)
end

# The SU(3) glueballs ladder system in the minimial truncation
if true
    # Important functions
    function ketbra(m::Integer, n::Integer, q::Integer)
        M = zeros(Int64, q, q)
        M[m, n] = 1
        return SparseMatrixCSC(M)
    end
    global lambda_from_g = g -> g^4 / (1 + g^4)
    function g_from_lambda(λ::Real)
        if λ < 0 || λ >= 1
            error("Invalid λ value")
        end
        return (λ/(1 - λ))^(1/4)
    end

    # Parameters
    global L0 = 3
    global Li = 20
    global d = 3 # The local dimension
    global λ = 0.1

    # Local operators
    global id = operator_identity(SparseMatrixCSC, 3)
    tau = SparseMatrixCSC([0 1 0; 0 0 1; 1 0 0])
    cc = SparseMatrixCSC([1 0 0; 0 0 1; 0 1 0])
    kb(n,m) = ketbra(n,m,3)
    casimir(n) = 2 // 3 * (3 - n) * n
    global E2two = spzeros(9, 9)
    for n1 in 0:2
        for n2 in 0:2
            i1 = n1 + 1
            i2 = n2 + 1
            m = mod(n1 - n2, 0:2)
            global E2two += (casimir(n1) + casimir(n2) + casimir(m)) * (kb(i1,i1) ⊗ kb(i2,i2))
        end
    end

    global E2 = 1/2 * (E2two ⊗ id) + 1/2 * (id ⊗ E2two) # Electric field (casimir) as a 3-local operator
    Dg1 = Diagonal([1, 1/√3, 1/3]) # First diagonal D terms for SU(3)
    Dg2 = Diagonal([1, 1/3, 1/√3]) # Second diagonal D terms for SU(3)
    Dg3 = Diagonal([1, 1/√3, 1/√3]) # Third diagonal D terms for SU(3)
    global Up = kron(Dg1, kb(3,1), Dg1) + kron(Dg2, kb(1,2), Dg2) + kron(Dg3, kb(2,3), Dg3)
    global H0 = λ * E2 - (1 - λ) * (Up + Up')

    # Operators of the pbc chain
    # global Hi = summation_local(H0, d, Li; pbc = true)
    # global Ti = operator_translation(SparseMatrixCSC, d, Li)
    # global Pi = operator_reflection(SparseMatrixCSC, d, Li)
end

#=
Hi_mat = summation_local(H0, d, Li; pbc = true)
Ti_mat = operator_translation(SparseMatrixCSC, d, Li)
disprel_mat = dispersion_relation(Hi_mat, Ti_mat, Li, nlevels = 8)
p1 = Plots.plot(disprel_mat)
=#

H0_mpo = mpo_from_matrix(Matrix(H0), d)
Hi_mpo = summation_local(H0_mpo, Li, pbc = true)
sitesl = siteinds_main(Hi_mpo)
Ti_mpo = operator_translation(MPO, d, Li)
replace_siteinds!(Ti_mpo, sitesl)
println("Test = ", norm(product(Hi_mpo, Ti_mpo) - product(Ti_mpo, Hi_mpo)))
disprel_mpo = dispersion_relation(Hi_mpo, nlevels = 5)
ωbs = pop_groundstate!(disprel) # The vacuum as Bloch state object
ωmps = wavefunction(ωbs)
E0 = energy(ωbs)
Plots.plot(disprel_mpo)

#=

Lbig = 50 # The length of the big system
χmax = 100 # The maximum bond dimension
jc = Int((Li - 1)/2 + 1) # The central site

# We compute the dispersion relation and we save the plot
disprel = dispersion_relation(Hi, Ti, Li; nlevels = 2)
plotdisprel = Plots.plot(disprel)

# We pop the groundstate and we compute its mps form
ωbs = pop_groundstate!(disprel) # The vacuum as Bloch state object
E0 = energy(ωbs)

# The vacuum in vector form
ωvec = wavefunction(ωbs)

# The vacuum in MPS form
sitesl = siteinds(d, Li)
ωmps = mps_from_vector(ωvec, d)
replace_siteinds!(ωmps, sitesl)

# We plot the dispersion relation
Plots.plot(plotdisprel)

=#