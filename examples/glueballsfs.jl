using Scattensor
using ITensors, ITensorMPS
using LinearAlgebra
using SparseArrays
using Optim
using Plots
using PlotlyJS
using KrylovKit
using Logging
using Revise
plotlyjs()

# We define the Irrep
struct Irp
    lab::Int
end

function Base.string(irp::Irp)
    if irp.lab == 0
        return "•"
    elseif irp.lab == 1
        return "▼"
    elseif irp.lab == 2
        return "▲"
    end
end

# We define the conjugate of irrep
function Base.adjoint(a::Irp)
    if a.lab == 0
        return Irp(0)
    elseif a.lab == 1
        return Irp(2)
    elseif a.lab == 2
        return Irp(1)
    end
end

# We define the fusion between irreps
function Base.:*(a::Irp, b::Irp)
    return Irp(mod((a.lab + b.lab), 3))
end

function quadratic_casimir(a::Irp)
    if a.lab == 0 
        return 0
    else
        return 4/3
    end
end

# We define the junction
struct Junc
    leg1::Irp
    leg2::Irp
    leg3::Irp
end

function isgauginv(jc::Junc)
    return (jc.leg1 * jc.leg2 * jc.leg3) == Irp(0)
    
end

function corner(jc::Junc)
    coeff = 0.0
    if jc == Junc(Irp(0), Irp(0), Irp(0)) || jc == Junc(Irp(0), Irp(2), Irp(1)) || jc == Junc(Irp(0), Irp(1), Irp(2))
        coeff = 1.0
    elseif jc == Junc(Irp(2), Irp(1), Irp(0)) || jc == Junc(Irp(1), Irp(0), Irp(2))
        coeff = 1.0 / sqrt(3)
    else
        coeff = 1.0 / sqrt(sqrt(3))
    end
    if jc == Junc(Irp(0), Irp(2), Irp(1)) || jc == Junc(Irp(1), Irp(1), Irp(1))
        coeff *= 1
    end
    return Junc(jc.leg1, jc.leg2 * Irp(2), jc.leg3 * Irp(1)), coeff
end

struct Plaq
    junc1::Junc
    junc2::Junc
    junc3::Junc
    junc4::Junc
end

function Base.print(plaq::Plaq)
    println()
    j1l1 = plaq.junc1.leg1
    j1l2 = plaq.junc1.leg2
    j1l3 = plaq.junc1.leg3
    j2l1 = plaq.junc2.leg1
    j2l2 = plaq.junc2.leg2
    j2l3 = plaq.junc2.leg3
    j3l1 = plaq.junc3.leg1
    j3l2 = plaq.junc3.leg2
    j3l3 = plaq.junc3.leg3
    j4l1 = plaq.junc4.leg1
    j4l2 = plaq.junc4.leg2
    j4l3 = plaq.junc4.leg3
    println("$(j1l1)", " ", "$(j1l2)", "$(j4l3)", " ", "$(j4l1)")
    println(" ", "$(j1l3)", " ", " ", "$(j4l2)", " ")
    println(" ", "$(j2l2)", " ", " ", "$(j3l3)", " ")
    println("$(j2l1)", " ", "$(j2l3)", "$(j3l2)", " ", "$(j3l1)")
    println()
end

function plaquetteoperator(plaq::Plaq)
    junc1, coeff1 = corner(plaq.junc1)
    junc2, coeff2 = corner(plaq.junc2)
    junc3, coeff3 = corner(plaq.junc3)
    junc4, coeff4 = corner(plaq.junc4)
    newPlaq = Plaq(junc1, junc2, junc3, junc4)
    coeff = coeff1 * coeff2 * coeff3 * coeff4
    return newPlaq, coeff
end

function electric_energy(plaq::Plaq)
    term1 = quadratic_casimir(plaq.junc1.leg2)
    term2 = quadratic_casimir(plaq.junc2.leg3)
    term3 = 0.5 * quadratic_casimir(plaq.junc2.leg2)
    term4 = 0.5 * quadratic_casimir(plaq.junc3.leg3)
    return term1 + term2 + term3 + term4
end

# We compute the basis of plaquettes
basis = []
for ia in 0:2
    for ib in 0:2
        for ic in 0:2
            la = Irp(ia)
            lb = Irp(ib)
            lc = Irp(ic)
            junc1 = Junc(la', lb, la * lb')
            junc2 = Junc(la, la' * lb, lb')
            junc3 = Junc(lc', lb, lb' * lc)
            junc4 = Junc(lc, lb * lc', lb')
            # We test that is all gauge invariant
            if !isgauginv(junc1) || !isgauginv(junc2) || !isgauginv(junc3) || !isgauginv(junc4)
                error("Not gauge invariant junction")
            end
            plaq = Plaq(junc1, junc2, junc3, junc4)
            push!(basis, plaq)
        end
    end
end

N = length(basis)
Uplaq = spzeros(N, N)
E2 = spzeros(N, N)
for ix in 1:N
    E2[ix, ix] = electric_energy(basis[ix])
    for iy in 1:N
        plaq, coeff = plaquetteoperator(basis[ix])
        if plaq == basis[iy]
            Uplaq[ix, iy] = coeff
        end
    end
end

HB = Uplaq + Uplaq'

function format_matrix_for_mathematica(matrix::AbstractMatrix)
    formatted_rows = ["{" * join(row, ", ") * "}" for row in eachrow(matrix)]
    return "{" * join(formatted_rows, ", ") * "}"
end

function matrix_to_mpo(A::Matrix, d::Int, N::Int; cutoff=1e-8)
    @assert size(A, 1) == d^N && size(A, 2) == d^N "Matrix dimensions must be d^N x d^N"
    sites = [Index(d, "Site, n=$i") for i in 1:N]  # For local dimension d=2, use "Qubit"; for d=3, use appropriate site type
    T = reshape(A, (repeat([d], 2N)...))
    IT = ITensor(T, (sites..., prime.(sites)...))
    W = MPO(IT, sites; cutoff=cutoff)
    return W, sites
end

function vector_to_mps(v::Vector, d::Int, N::Int; cutoff=1e-8, maxdim=10)
    @assert length(v) == d^N "Vector length must be d^N"
    sites = [Index(d, "Site, n=$i") for i in 1:N]
    T = reshape(v, (repeat([d], N)...))
    IT = ITensor(T, sites...)
    psi = MPS(IT, sites; cutoff=cutoff, maxdim=maxdim)
    return psi, sites
end

"""
Constructs the reflection operator matrix for a 1D system with L sites and local dimension d.
"""
function reflection_operator(L::Int, d::Int)
    # Generate all possible basis states as integer arrays
    basis_states = collect(Iterators.product(ntuple(_ -> 0:d-1, L)...))
    dim_H = d^L  # Total Hilbert space dimension
    R = spzeros(Int, dim_H, dim_H)  # Initialize reflection matrix
    
    # Map basis states to their reflected counterparts
    state_to_index = Dict(state => i for (i, state) in enumerate(basis_states))
    
    for (i, state) in enumerate(basis_states)
        reflected_state = reverse(state)  # Reflect the state
        j = state_to_index[reflected_state]  # Get index of the reflected state
        R[j, i] = 1  # Set matrix element
    end
    
    return R
end

"""
Constructs the reflection operator matrix for a 1D system with L sites and local dimension d.
"""
function translation_operator(L::Int, d::Int)
    # Generate all possible basis states as integer arrays
    basis_states = collect(Iterators.product(ntuple(_ -> 0:d-1, L)...))
    dim_H = d^L  # Total Hilbert space dimension
    T = spzeros(Int, dim_H, dim_H)  # Initialize reflection matrix
    
    # Map basis states to their reflected counterparts
    state_to_index = Dict(state => i for (i, state) in enumerate(basis_states))
    
    for (i, state) in enumerate(basis_states)
        reflected_state = Tuple(circshift([i for i in state], 1))  # Reflect the state
        j = state_to_index[reflected_state]  # Get index of the reflected state
        T[j, i] = 1  # Set matrix element
    end
    
    return T
end

function glueball_spectrum_lanczos(L, λ)
    gE = λ
    gB = 1 - λ
    eps = (1 - λ)/λ

    # We define the local Hamiltonian
    H0 = LocalOperator(gE * E2 - gB * (Uplaq + Uplaq'), [3, 3, 3], "h")
    drel = disprel(H0, L, nlevels = 30)

    println("Grounstate energy = ", energy(getgroundstate(drel)))

    Plots.plot(drel)

    # Groundstate energy
    # We take only the k = 0
    # for el in drel
    #     if el.kfraction == 0
    #         println((energy(el) - 16/3 - eps)/(0.01)^2, ", ", (energy(el) - 16/3 + eps)/(0.01)^2)
    #     end
    # end

    return drel
end


# Define a system with qutrits
ITensors.space(::SiteType"glueball") = 3
bv1 = [1 0 0]
bv2 = [0 1 0]
bv3 = [0 0 1]
operatorc = [0  1  0
             0  0  1
             1  0  0]
operatorn0 =[1  0  0
             0  0  0
             0  0  0]
operatorn1 =[0  0  0
             0  1  0
             0  0  0]
operatorn2 =[0  0  0
             0  0  0
             0  0  1]
ITensors.op(::OpName"cdag", ::SiteType"glueball") = operatorc
ITensors.op(::OpName"c", ::SiteType"glueball") = operatorc'
ITensors.op(::OpName"c+", ::SiteType"glueball") = (operatorc + operatorc') / 2
ITensors.op(::OpName"c-", ::SiteType"glueball") = (operatorc - operatorc') / 2
ITensors.op(::OpName"n0", ::SiteType"glueball") = operatorn0
ITensors.op(::OpName"n1", ::SiteType"glueball") = operatorn1
ITensors.op(::OpName"n2", ::SiteType"glueball") = operatorn2
ITensors.op(::OpName"n+", ::SiteType"glueball") = (operatorn1 + operatorn2) / 2
ITensors.op(::OpName"n-", ::SiteType"glueball") = (operatorn1 - operatorn2) / 2
ITensors.op(::OpName"bas11", ::SiteType"glueball") = bv1'*bv1
ITensors.op(::OpName"bas12", ::SiteType"glueball") = bv1'*bv2
ITensors.op(::OpName"bas13", ::SiteType"glueball") = bv1'*bv3
ITensors.op(::OpName"bas21", ::SiteType"glueball") = bv2'*bv1
ITensors.op(::OpName"bas22", ::SiteType"glueball") = bv2'*bv2
ITensors.op(::OpName"bas23", ::SiteType"glueball") = bv2'*bv3
ITensors.op(::OpName"bas31", ::SiteType"glueball") = bv3'*bv1
ITensors.op(::OpName"bas32", ::SiteType"glueball") = bv3'*bv2
ITensors.op(::OpName"bas33", ::SiteType"glueball") = bv3'*bv3

function ↻(n::Integer, m::Integer)::Integer
    n > 0 ? (n-1)%m + 1 : m + n%m
end

# GLOBAL VARIABLES
Lglobal = 20
λglobal = 0.3
sitesglobal = siteinds("glueball", Lglobal)

function hamiltonian_mpo(; sites = sitesglobal, λ = λglobal)
    L = length(sites)
    ampo = AutoMPO()
    H0r = reshape(λ * E2 - (1 - λ) * (Uplaq + Uplaq'), (3,3,3,3,3,3))
    for j in 1:L
        ranges = [1:3 for _ in 1:6]
        for tuple in Iterators.product(ranges...)
            H0rvalue = H0r[tuple[1], tuple[2], tuple[3], tuple[4], tuple[5], tuple[6]]
            if H0rvalue != 0
                oprj = "bas$(tuple[1])$(tuple[4])"
                oprjp1 = "bas$(tuple[2])$(tuple[5])"
                oprjp2 = "bas$(tuple[3])$(tuple[6])"
                add!(ampo, H0rvalue, oprj, j↻L , oprjp1, (j+1)↻L, oprjp2, (j+2)↻L)
            end
        end
    end
    return MPO(ampo, sites)
end

function local_hamiltonian_mpo(j; sites = sitesglobal, λ = λglobal)
    L = length(sites)
    ampo = AutoMPO()
    H0r = reshape(λ * E2 - (1 - λ) * (Uplaq + Uplaq'), (3,3,3,3,3,3))
    ranges = [1:3 for _ in 1:6]  # Create a list of ranges
    for tuple in Iterators.product(ranges...)
        H0rvalue = H0r[tuple[1], tuple[2], tuple[3], tuple[4], tuple[5], tuple[6]]
        if H0rvalue != 0
            oprj = "bas$(tuple[1])$(tuple[4])"
            oprjp1 = "bas$(tuple[2])$(tuple[5])"
            oprjp2 = "bas$(tuple[3])$(tuple[6])"
            add!(ampo, H0rvalue, oprj, j↻L , oprjp1, (j+1)↻L, oprjp2, (j+2)↻L)
        end
    end
    return MPO(ampo, sites)
end

Hmpo = hamiltonian_mpo()
Hlocmpo = [local_hamiltonian_mpo(j) for j in 1:Lglobal]

Eng0, ψ0 = dmrg(Hmpo,
                randomMPS(sitesglobal; linkdims=100);
                nsweeps=10, 
                maxdim=[100, 200], 
                cutoff=1e-12)

println("Ground state energy = $Eng0")

function apply_localop_to_mps(localopname, mps, j)
    newmps = deepcopy(mps)
    s = siteind(newmps, j)
    Opt = op(localopname, s)
    new_tensor = Opt * newmps[j]
    noprime!(new_tensor)
    newmps[j] = new_tensor
    return newmps
end

function wavepacketcreator(jbar, kbar, sigma; sites = sitesglobal)
    L = length(sites)
    autmpo = AutoMPO()
    for j in 1:L  # Ensure we stay within bounds
        coeff = exp(im * kbar * j) * exp(-(j - jbar)^2 / (2 * sigma))
        if abs(coeff) > 10^(-15)
            add!(autmpo, coeff, "c+", j↻L)  # Example spin interaction
        end
    end
    return MPO(autmpo, sites)
end

function energyprofile(ψmps, ψzero; L = Lglobal)
    return [real(inner(ψmps', Hlocmpo[j], ψmps) - inner(ψzero', Hlocmpo[j], ψzero)) for j in 1:(L-2)]
end

function particleprofile(ψmps; L = Lglobal)
    ret = []
    for j in 2:(L-1)
        ψpart = apply_localop_to_mps("n+", ψmps, j)
        push!(ret, abs(inner(ψmps', ψpart)))
    end
    return ret
end

ψ1 = apply_localop_to_mps("c+", ψ0, 10)
ψ1 = normalize(ψ1)
# print("1-particle energy = ", inner(ψ1', H, ψ1))

# ψ1 = normalize(apply(wavepacketcreator(10, π/2, 1.5), ψ0))
# Plots.plot(energyprofile(ψ1, ψ0))

function tdvp_time_evolution(H::MPO, ψ::MPS, ψ₀::MPS, dt::Float64, Δt::Float64; L = Lglobal)

    N = Integer(round(Δt / dt)) # Number of steps
    # We create the list of times 
    times = [n * dt for n in 1:N]
    positions = [x for x in 2:L-1]

    # We create a 2D array E of dimensions (N, L-2)
    En = zeros(N, L-2)
    Pp = zeros(N, L-2)
    EnLog = zeros(N, L-2)
    Linkdims = zeros(N, L-2)
    Maxlinkdim = zeros(N)

    for n in 1:N
        println("Step ", n, " of ", N)

        timeElapsed = @elapsed begin
            ψ = tdvp(H, -im * dt, ψ, maxdim = 10)
            normalize!(ψ)
            # print("maxdim1 = ", maxlinkdim(ψ))
            # orthogonalize!(ψ, 1, maxdim = 20)
            # orthogonalize!(ψ, L, maxdim = 20)
            print("maxdim = ", maxlinkdim(ψ))
            Eprofile = energyprofile(ψ, ψ₀)
            Pprofile = particleprofile(ψ)
            for j in 2:(L-1)
                En[n,j-1] = Eprofile[j-1]
                Pp[n,j-1] = Pprofile[j-1]
                EnLog[n,j-1] = log10(abs(Eprofile[j-1]))
                Linkdims[n,j-1] = linkdim(ψ,j)
                Maxlinkdim[n] = maxlinkdim(ψ)
            end
        end

        # Print information
        println("Step ", n, " of ", N, " completed.")
        println("Maximum link dimension of MPS: ", maxlinkdim(ψ))
        println("Relative Energy: ", real(inner(ψ', H, ψ) - inner(ψ₀', H, ψ₀)))
        println("Estimated time remained: ", timeElapsed * (N - n), " seconds.")
    end

    # heatmap(positions, times, En)
    # savefig("simulation_output/En.png")
    # heatmap(positions, times, EnLog)
    # savefig("simulation_output/EnLog.png")
    # heatmap(positions, times, Linkdims)
    # savefig("simulation_output/Linkdims.png")
    # titleString = "L = $L, λ = $λ, jⁱᵐᵖ = $jⁱᵐᵖ, jᵂᴾ = $jᵂᴾ, kᵂᴾ = $kᵂᴾ,
                    # σᵂᴾ = $σᵂᴾ, ϵᵂᴾ = $ϵᵂᴾ, Nᴰᴹᴿᴳ = $Nᴰᴹᴿᴳ, χᴰᴹᴿᴳ ≤ $χᴰᴹᴿᴳ,
                    # ϵᴰᴹᴿᴳ = $ϵᴰᴹᴿᴳ, χᵀᴰⱽᴾ ≤ $χᵀᴰⱽᴾ, ϵᵀᴰⱽᴾ = $ϵᵀᴰⱽᴾ, dt = $dt, Δt = $Δt"
    # titleString = "Graph"
    # l = @layout [grid(2,2); b{0.2h}]
    # title = Plots.plot(title = titleString, grid = false, showaxis = false)
    # Plots.plot(p1, p2, p3, p4, title, layout = l, dpi = 300, size=(1000, 1000))

    return Dict(
        "E" => Plots.heatmap(positions, times, En),
        "Pprofile" => Plots.heatmap(positions, times, Pp),
        "LogE" => Plots.heatmap(positions, times, EnLog), 
        "Linkdims" => Plots.heatmap(positions, times), 
        "Maxlinkdims" => Plots.plot(times, Maxlinkdim))
end

myplots = tdvp_time_evolution(Hmpo, ψ1, ψ0, 0.1, 10.0)