using Makie
using MixedModels
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

m1 = fit(
    MixedModel,
    @formula(1000/reaction ~ 1 + days + (1+days|subj)),
    MixedModels.dataset(:sleepstudy),
)

@testset "ranefinfo" begin
    reinfo = ranefinfo(m1)
    @test isone(length(reinfo))
    @test keys(reinfo) == (:subj, )
    re1 = first(reinfo)
    @test isa(re1, RanefInfo)
    @test re1.cnames == ["(Intercept)", "days"]
    @test first(re1.levels) == "S308"
    @test first(re1.ranef) ≈ -0.081830283935302 atol=1.0e-5
    @test first(re1.stddev) ≈ 0.14035486016644683 atol=1.0e-5
    f = Figure();
    @test f.content == Any[]
    caterpillar!(f, re1);
    @test length(f.content) == 2
    @test isa(first(f.content), Axis)
    @test isone(f.layout.nrows)
    @test f.layout.ncols == 2
end

@testset "shrinkageplot" begin
    f = shrinkageplot(m1)
    @test isone(length(f.content))
end

@testset "utilities" begin
    ppts = MixedModelsMakie.ppoints(64)
    @test length(ppts) == 64
    @test first(ppts) ≈ inv(128)
end
