"""
Interpolates a state 𝛙₁ from a state 𝛙₀, using a linear combinations of compositions
of local operators 𝐋ⱼ.
"""
function interpolate(
    𝛙₁::Vector{ComplexF64}, # The vector to interpolate
    𝛙₀::Vector{ComplexF64}, # The vector from which 𝛙₁ is interpolated
    𝐓::Matrix{ComplexF64}, # The translation operator of the system
    j₁::Int64, # The position of the left border of the interpolation range
    j₂::Int64, # The position of the right border of the interpolation range
    𝐋::Vector{Matrix{ComplexF64}}, # A list of d local operators
    d::Int64, # The dimension of the local space
    L::Int64, # The number of sites of the chain
)
    # We define the identity
    Ide = Matrix{ComplexF64}(I, d, d)

    function generate_𝐋₁(n)
        𝐀 = 𝐋[n]
        for _ in 2:L
            𝐀 = 𝐀 ⊗ Ide
        end
        return 𝐀
    end

    𝐋₁ = [generate_𝐋₁(n) for n in 1:d]

    # Apply 𝐋ⱼ to the state 𝛙 in the most efficient way (translating the state instead 
    # of the operator)
    function apply(n, j, 𝛙)
        𝛟 = deepcopy(𝛙)
        for _ in 1:j
            𝛟 = 𝐓' * 𝛟
        end
        𝛟 = 𝐋₁[n+1] * 𝛟
        for _ in 1:j
            𝛟 = 𝐓 * 𝛟
        end
        return 𝛟
    end

    function apply_to_array(n, j, arr)
        arr′ = deepcopy(arr)
        arr′[j] = n
        return arr′
    end

    𝚿 = [𝛙₀]
    A = [zeros(L)]

    # Find all the possible combinations of the operators starting from 𝛙₀
    # We also find an array representing these applications
    for j in j₁:j₂
        𝚿 = [apply(n,j,𝛙) for n in 0:(d-1) for 𝛙 in 𝚿]
        A = [apply_to_array(n,j,a) for n in 0:(d-1) for a in A]
    end

    # We compute the Hessian matrix
    Hs = [𝚿[α]' * 𝚿[β] for α in eachindex(𝚿), β in eachindex(𝚿)]
    b = [𝚿[α]' * 𝛙₁ for α in eachindex(𝚿)]
    println("Dimension of the Hessian matrix: ", size(Hs))
    println("Rank of the Hessian matrix: ", rank(Hs))

    # We compute the coefficients, calculating the pseudo-inverse of 𝚿
    println("Computing the inverse of the Hessian matrix...")
    c = pinv(Hs) * b

    # We compute the infidelity
    iFi = 1 + sum(c[α]' * c[β] * Hs[α,β] for α in eachindex(c), β in eachindex(c))
                 - sum((c[α]' * b[α] + c[α] * b[α]') for α in eachindex(c))
    println("Infidelity: ", iFi)

    # We create a list of the tuples (coefficients, A)
    listt = [(c[α], A[α]) for α in eachindex(c)]

    # We sort the list by the coefficients
    sort!(listt, by = x->abs(x[1]), rev = true)

    # We print the first 20 elements of listt in a nice way
    # println()
    # println("The first 20 elements are:")
    # println()
    # for α in 1:20
        # println(listt[α][2], " ", abs(listt[α][1]))
    # end

    # we create the array of the abs of the coefficients sorted
    xax = [i for i in 1:length(listt)]
    cp = [abs(x[1]) for x in listt]

    plot(xax, cp, xaxis=:log, yaxis=:log)
    savefig("data/coefficients.png")
end
export interpolate