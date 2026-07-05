# Arctan profile fit: recovery of synthetic resonance parameters.
@testset "fits" begin
    @test wrap_mod_pi(3.0π + 0.2) ≈ 0.2 atol = 1e-12
    @test wrap_mod_pi(-0.4) ≈ -0.4 atol = 1e-12

    x0, w, A = 2.84, 0.05, 1.0
    xs = collect(range(2.5, 3.3; length = 25))
    ys = 0.3 .- 0.8 .* (xs .- 2.9) .+ A .* atan.(2 .* (xs .- x0) ./ w)
    fit = fit_arctan(xs, ys; bg_order = 1)
    @test fit.x0 ≈ x0 atol = 0.01
    @test fit.w ≈ w rtol = 0.3
    @test fit.A ≈ A rtol = 0.05
    @test fit.rms < 0.02
end
