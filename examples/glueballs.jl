
using LinearAlgebra
using SparseArrays
using Plots

using Scattensor
using Revise

function H0_glueballs(; gE = 0, gB = 0)
    
    # link states: 0 trivial, 1 fundamental, 2 antifundamental
    # Positive direction: exiting from junction
    # 0 trivial, 1 fundamental, 2 antifundamental
    # A 3-junction is a tuple.
    # The constructor of a three-junction
    # can represent a top T junction or a bottom T junction
    function junction(l1::Int, l2::Int, l3::Int)
        return (mod(l1, 3), mod(l2, 3), mod(l3, 3))
    end

    function isgauginv(j)
        return (mod(j[1] + j[2] + j[3], 3) == 0)
    end

    # va = 1.0 / sqrt(3)
    # vb = 1.0 / sqrt(sqrt(27))
    # vc = 1.0 / 3.0

    # Coefficients of UR in the T junction (of basis by Pietro)
    # C1 = va
    # C2 = va
    # C3 = va
    # C4 = -vb
    # C5 = vc
    # C6 = -vb
    # C7 = vb
    # C8 = vc
    # C9 = vb

    C1 = 1
    C2 = 1
    C3 = 1
    C4 = 1
    C5 = 1/3
    C6 = 1
    C7 = 1
    C8 = 1/3
    C9 = 1

    # U right for gauge invariant top-T-junction (Left, Right, Bottom links)
    function URT(j::Tuple)
        j′ = junction(j[1] + 1, j[2], j[3] - 1)
        if j == (0,0,0)
            coeff = C1
        elseif j == (0,2,1)
            coeff = C2
        elseif j == (0,1,2)
            coeff = C3
        elseif j == (1,1,1)
            coeff = C4
        elseif j == (1,0,2)
            coeff = C5
        elseif j == (1,2,0)
            coeff = C6
        elseif j == (2,2,2)
            coeff = C7
        elseif j == (2,1,0)
            coeff = C8
        elseif j == (2,0,1)
            coeff = C9
        end
        return (j′, coeff)
    end

    # U left for gauge invariant top-T-junction (Left, Right, Bottom links)
    function ULT(j::Tuple)
        j′ = junction(j[1], j[2] - 1, j[3] + 1)
        if j == (0,0,0)
            coeff = C1
        elseif j == (0,2,1)
            coeff = C8
        elseif j == (0,1,2)
            coeff = C6
        elseif j == (1,1,1)
            coeff = C4
        elseif j == (1,0,2)
            coeff = C2
        elseif j == (1,2,0)
            coeff = C9
        elseif j == (2,2,2)
            coeff = C7
        elseif j == (2,1,0)
            coeff = C5
        elseif j == (2,0,1)
            coeff = C3
        end
        return (j′, coeff)
    end

    # U left for gauge invariant bottom-T-junction (Left, Right, Top links)
    function ULB(j::Tuple)
        j′ = junction(j[1] - 1, j[2], j[3] + 1)
        if j == (0,0,0)
            coeff = C1
        elseif j == (0,2,1)
            coeff = C9
        elseif j == (0,1,2)
            coeff = C5
        elseif j == (1,1,1)
            coeff = C4
        elseif j == (1,0,2)
            coeff = C3
        elseif j == (1,2,0)
            coeff = C8
        elseif j == (2,2,2)
            coeff = C7
        elseif j == (2,1,0)
            coeff = C6
        elseif j == (2,0,1)
            coeff = C2
        end
        return (j′, coeff)
    end

    # U right for gauge invariant bottom-T-junction (Left, Right, Top links)
    function URB(j::Tuple)
        j′ = junction(j[1], j[2] + 1, j[3] - 1)
        if j == (0,0,0)
            coeff = C1
        elseif j == (0,2,1)
            coeff = C3
        elseif j == (0,1,2)
            coeff = C2
        elseif j == (1,1,1)
            coeff = C4
        elseif j == (1,0,2)
            coeff = C6
        elseif j == (1,2,0)
            coeff = C5
        elseif j == (2,2,2)
            coeff = C7
        elseif j == (2,1,0)
            coeff = C9
        elseif j == (2,0,1)
            coeff = C8
        end
        return (j′, coeff)
    end

    # A plaquette is a (named) tuple of four junctions:
    # (TL = (top-left), TR = (top-right), BL = (bottom-left), BR = (bottom-right))

    function U(p::NamedTuple)
        p′ = (TL = URT(p[:TL])[1], TR = ULT(p[:TR])[1], BL = URB(p[:BL])[1], BR = ULB(p[:BR])[1])
        coeff = URT(p[:TL])[2] * ULT(p[:TR])[2] * URB(p[:BL])[2] * ULB(p[:BR])[2]
        return (p′, coeff)
    end

    function longpol(p::NamedTuple)
        return mod((p[:TL][2] + p[:BL][2]),3)
    end

    function isconsistent(p::NamedTuple)
        c1 = (mod(p[:TL][2] + p[:TR][1], 3) == 0)
        c2 = (mod(p[:TL][3] + p[:BL][3], 3) == 0)
        c3 = (mod(p[:BR][1] + p[:BL][2], 3) == 0)
        c4 = (mod(p[:BR][3] + p[:TR][3], 3) == 0)
        return c1 && c2 && c3 && c4
    end

    # We construct the basis of the space of the top-junctions
    juncset = []
    for l in [0,1,2] 
        for r in [0,1,2]
            for v in [0,1,2]
                j = junction(l, r, v)
                if isgauginv(j)
                    push!(juncset, j)
                end
            end
        end
    end

    # We construct the basis of the space of plaquettes
    basis = []
    for tl in juncset
        for tr in juncset
            for bl in juncset
                for br in juncset
                    p = (TL = tl, TR = tr, BL = bl, BR = br)
                    if isconsistent(p)
                        if longpol(p) == 0
                            push!(basis, p)
                        end
                    end
                end
            end
        end
    end

    # We define the dimension of the space of plaquettes
    dim = length(basis)

    # We define a function that assigns a unique number to each plaquette
    function chainpos(p::NamedTuple)
        v1 = p[:BL][1]
        v2 = p[:BL][2]
        v3 = p[:BR][2]
        return 3^2 * v1 + 3^1 * v2 + 3^0 * v3 + 1
    end

    # We sort the basis by chainpos
    sort!(basis, by = p -> chainpos(p))

    # We constructq the matrix of the operator U
    Uplaq = spzeros(dim, dim)
    for p in basis
        Uplaq[chainpos(p), chainpos(U(p)[1])] = U(p)[2]
    end

    # We plot the matrix Uplaquette
    # Plots.heatmap(Matrix(Uplaq))
    # Plots.savefig("Uplaq.png")

    # We plot the matri Uplaquette * Uplaquette'
    # Plots.heatmap(Matrix(Uplaq * Uplaq'))
    # Plots.savefig("UplaqUplaq.png")


    Casimir2 = 4/3
    E2 = spzeros(dim, dim)
    for p in basis
        function f(x)
            if x == 0
                return 0
            elseif x == 1
                return Casimir2
            elseif x == 2
                return Casimir2
            else
                error("Invalid value")
            end
        end
        electric = f(p[:TL][2]) + f(p[:BL][2]) + 1/2 * f(p[:TL][3]) + 1/2 * f(p[:TR][3])
        E2[chainpos(p), chainpos(p)] = electric
    end

    return gE * E2 - gB * (Uplaq + Uplaq')
end

# We define the local Hamiltonian
λ = 0.0
H0matrix = H0_glueballs(gE = λ, gB = 1 - λ)

H0 = LocalOperator(H0matrix, [3, 3, 3], "h")

println("Operator H0 created")

# drel = disprel(H0, 7)