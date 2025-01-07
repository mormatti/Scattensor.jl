using LinearAlgebra
using .Utils.BasicTools
using .Utils.LinearAlgebraExtensions

# Method from analytical formulas. Not use this if possible.
# Use just numerical minimization.

n = 5

if true
    # We create a random 5x5 hermitian matrix
    Λ = 3 * ((rand(Float64, n, n) - 0 * ones(Float64, n, n)) + im * (rand(Float64, n, n) - 0 * ones(Float64, n, n)))
    Λ = (Λ + Λ')/2
    display(Λ)
end

let 
    z = rand(ComplexF64, n)
    z = z ./ abs.(z)

    # We iterate a steepest descent
    for _ in 1:3
        for γ in eachindex(z)
            s = sum([Λ[α,γ] * z[α] for α in eachindex(z)]) - z[γ] * Λ[γ,γ]
            z1γ = sqrt(s / s')
            z2γ = -z1γ
            z1 = deepcopy(z)
            z2 = deepcopy(z)
            z1[γ] = z1γ
            z2[γ] = z2γ
            z = (abs(z1' * Λ * z1) > abs(z2' * Λ * z2)) ? z1 : z2
            println("f = ", z' * Λ * z)
        end
    end
    
    z = z / z[1]
    println("Output:")
    # Ee write the first 5 elements of z 
    println((z[1:5])')
end