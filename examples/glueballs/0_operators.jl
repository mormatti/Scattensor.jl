module Glueballs

    using ITensorMPS
    using LinearAlgebra
    using SparseArrays
    using Scattensor

    # The structure for the Irrep
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

    basislengthconf = length(basis)
    Uplaq = spzeros(basislengthconf, basislengthconf)
    E2 = spzeros(basislengthconf, basislengthconf)
    for ix in 1:basislengthconf
        E2[ix, ix] = electric_energy(basis[ix])
        for iy in 1:basislengthconf
            plaq, coeff = plaquetteoperator(basis[ix])
            if plaq == basis[iy]
                Uplaq[ix, iy] = coeff
            end
        end
    end

    """ The 3-local plaquette operator in matrix form (without pre-factor).
    """
    function plaquette_operator()
        return Uplaq
    end
    export plaquette_operator

    """ The 3-local operator H_B in matrix form (without pre-factor).
    """
    function magnetic_operator()
        return Uplaq + Uplaq'
    end
    export magnetic_operator

    """ The 3-local operator H_E in matrix form (without pre-factor).
    """
    function electric_operator()
        return E2
    end
    export electric_operator

    """ The 1-local operator n which couts the state n (1 = singlet, 
    2 = fundamental, 3 = antifundamental)
    """
    function charge_operator(n::IntType) where {IntType <: Integer}
        op = zeros(3,3)
        op[1,1] = 0
        op[2,2] = 1
        op[3,3] = -1
        return op
    end
    export charge_operator

    """ The 1-local operator b which increases a state but not cyclically.
    """
    function increase_operator()
        op = zeros(3,3)
        op[1,2] = 1
        op[2,3] = 1
        return op
    end
    export increase_operator

    """ The 1-local operator which increases a state cyclically.
    """
    function cyclic_operator()
        op = zeros(3,3)
        op[1,3] = 1
        op[2,3] = 1
        op[3,1] = 1
        return op
    end
    export cyclic_operator

    """ The local operator c which inverts (conjugates) the representation.
    """
    function conjugation_operator()
        op = zeros(3,3)
        op[1,1] = 1
        op[2,3] = 1
        op[3,2] = 1
        return op
    end
    export conjugation_operator

    """ Get an element of the basis in bra form: [0 0 1], [0 1 0], [0 0 1].
    """
    function bra(n::IntType) where {IntType <: Integer}
        if n == 1
            return [1 0 0]
        elseif n == 2
            return [0 1 0]
        elseif n == 3
            return [0 0 1]
        else
            error("Invalid number n for the bra")
        end
    end
    export bra

    """ The 3-local projector for the particle 0++
    """
    function projectr(v::Vector{Int64})
        # We check that all entries are between 1 and 3
        for vel in v
            if vel > 3 || vel < 1
                error("Entries must be between 1 and 3")
            end
        end

        vecof3 = [3 for _ in eachindex(v)]
        state = zeros(vecof3...)
        state[v...] = 1
        reshaped = reshape(state, prod(vecof3))
        proj = kron(reshaped', reshaped)
        mps = convert_to(MPO, proj, 3, length(v); cutoff = 1e-12)
        return mps
    end
    export projectr
end