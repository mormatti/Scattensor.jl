using Scattensor

@testset "Scattensor.jl" begin
    println("Testing Scattensor.jl")

    T = translation_operator(3, 2)
    println(T)
end