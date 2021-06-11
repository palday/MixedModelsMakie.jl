using MixedModelsMakie
using MixedModelsMakie.MixedModels
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
    @test a ≈ first(result) atol=0.05
    @test b ≈ last(result) atol=0.05
end

@testset "Shrinkage" begin
    
    m1 = let form = @formula(1000/reaction ~ 1 + days + (1+days|subj))
        fit(MixedModel, form, MixedModels.dataset(:sleepstudy))
    end
    shrk1 = shrinkage(m1)
    @test isone(length(shrk1))
    @test first(keys(shrk1)) == :subj
    glob = shrk1.subj.globalest
    @test isa(glob, NamedTuple)
    @test length(glob) == 2
    @test keys(glob) == (Symbol("(Intercept)"), :days)
    mxd = shrk1.subj.withinmxdtbl
    @test isa(mxd, NamedTuple)
    @test keys(mxd) == (:level, :mixed, :within)
end
