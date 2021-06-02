using MixedModelsMakie
using Random # we don't depend on exact PRNG vals, so no need for StableRNGs
using Test

@testset "There are no graphical tests" begin
    @test true
end

@testset "Simple linear regression" begin
    a, b = 1, 2
    x = collect(1:10)
    y = randn(MersenneTwister(42), 10) * 0.1
    @. y += a + b * x
    result = simplelinreg(x, y)
    @test result isa Tuple
    @test a ≈ result[1] atol=0.05
    @test b ≈ result[2] atol=0.05
end
