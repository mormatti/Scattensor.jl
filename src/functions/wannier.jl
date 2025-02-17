
function wannier(
    ψ::Vector{BlochState{Vector{ComplexType}}}, # The states of the band,
    H0::LocalOperator{MatrixType}, # The local hamiltonian of the system
    E0::Float64;
    ) where {
        ComplexType <: Complex,
        MatrixType <: Union{Hermitian, SparseMatrixCSC}
        }

    # This condition, in our case, is always satisfied
    L = length(ψ)

    # We first check that ψ is really a single band, i.e. that two different states
    # has different momentum.
    for i in 1:length(ψ)
        for j in 1:length(ψ)
            if i != j
                if ψ[i].kfraction == ψ[j].kfraction
                    error("ψ must be a single band.")
                end
            end
        end
    end

    # We define the constants
    H = Diagonal([energy(ψ[i]) for i in 1:length(ψ)])

    H -= E0 * I # We perform the shift

    T = Diagonal([exp(im * momentum(ψ[i])) for i in 1:length(ψ)])

    if !hasuniformlocaldims(H0)
        error("The local dimensions of the sites must be the same.")
    end
    d = H0.localdims[1]
    L0 = length(H0)
    H0ext = kron(H0.repr, generateidentity(MatrixType, d^(L-L0)))

    H1 = Hermitian([wavefunction(ψ[i1])' * H0ext * wavefunction(ψ[i2]) for i1 in 1:length(ψ), i2 in 1:length(ψ)])
    H1 -= E0 / L * I

    # We define the distance function
    # If L0 is even, we increment by 1/2 jc
    jc = Int64(floor((L+1)/2))
    if L0 % 2 == 0
        jc += 1/2
    end
    χ(j) = sin(π * (j - jc) / L)

    # We compute all the Hj, trH, X1 and X2
    Hj = [Hermitian(T^(-j+1) * H1 * T^(j-1)) for j in 1:L]

    trH = tr(H)
    X1 = sum(Hj[j] * χ(j) for j in 1:L)
    X2 = sum(Hj[j] * χ(j)^2 for j in 1:L)

    function funct(θ)
        w = exp.(im * θ)
        return Float64((1/trH) * real(w' * X2 * w)) - Float64((1/trH)^2 * real(w' * X1 * w)^2)
    end

    println()
    opt = optimize(funct, rand(Float64, length(ψ)), SimulatedAnnealing())
    opt = optimize(funct, Optim.minimizer(opt), NelderMead())

    θopt = Optim.minimizer(opt)
    wopt = exp.(im * θopt)
    engs = [real(wopt' * Hj[j] * wopt) for j in 1:L]

    Plots.plot(engs)

    return engs, opt

    # The following if you need the Hessian

    # function expval2(M, w)
    #     return expval(M, w)^2
    # end

    # function gradexpval(M, w)
    #     result = (im * w)' .* (M * conj(w))
    #     result += conj(result)
    #     return result
    # end

    # function gradexpval2(M, w)
    #     return 2 * expval(M, w) * gradexpval(M, w)
    # end

    # function hessianexpval(M, w)
    #     result = M .* (w * w') + M * w * w' .* I
    #     result += conj(result)
    #     return result
    # end

    # function hessianexpval2(M, w)
    #     return 2 * expval(M, w) * hessianexpval(M, w) + 2 * gradexpval(M, w) * gradexpval(M, w)'
    # end

    # function gradient!(storage, θ)
    #     w = exp.(im * θ)
    #     storage = (1/TrH) * gradexpval(X2, w) - (1/TrH)^2 * gradexpval2(X1, w)
    # end

    # function hessian!(storage, θ)
    #     w = exp.(im * θ)
    #     storage = (1/TrH) * hessianexpval(X2, w) - (1/TrH)^2 * hessianexpval2(X1, w)
    # end
end

# Step 3: l'utente conferma che la Wannier va bene, e da questa calcola l'operatore di
# creazione della particella, lo applica al vuoto DMRG e calcola le funzioni di
# correlazione da cui si possono estrarre le simulazioni di scattering.