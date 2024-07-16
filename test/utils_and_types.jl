@testset "utilities" begin
    ppts = MixedModelsMakie.ppoints(64)
    @test length(ppts) == 64
    @test first(ppts) ≈ inv(128)

    @test MixedModelsMakie.zquantile(0.025) ≈ -1.96 atol = 0.005
    @test MixedModelsMakie.zquantile(0.975) ≈ 1.96 atol = 0.005
    @test MixedModelsMakie.zquantile(0.50) ≈ 0
end

@testset "Simple linear regression" begin
    a, b = 1, 2
    n = 100
    x = 1:n
    y = randn(MersenneTwister(42), n) * 0.1
    @. y += a + b * x
    result = simplelinreg(x, y)
    @test result isa Tuple
    @test a ≈ result[1] atol = 0.05
    @test b ≈ result[2] atol = 0.05
end

m1 = fit(MixedModel,
         @formula(1000 / reaction ~ 1 + days + (1 + days | subj)),
         MixedModels.dataset(:sleepstudy); progress)

@testset "confint_table" begin
    wald = confint_table(m1, 0.68)
    bsamp = parametricbootstrap(MersenneTwister(42), 1000, m1; progress)
    boot = confint_table(bsamp, 0.68)

    @test wald.coefname == boot.coefname
    @test wald.estimate ≈ boot.estimate rtol = 0.05
    @test wald.lower ≈ boot.lower rtol = 0.05
    @test wald.upper ≈ boot.upper rtol = 0.05

    @test all(splat(isapprox),
              zip(MixedModelsMakie.confint_table(mr).estimate, fixef(mr)))

    @test fixefnames(mr) == MixedModelsMakie.confint_table(mr).coefname
    @test fixefnames(mr) == MixedModelsMakie.confint_table(br).coefname
end

@testset "ranefinfo" begin
    reinfo = ranefinfo(m1)
    @test isone(length(reinfo))
    @test keys(reinfo) == (:subj,)
    re1 = only(reinfo)
    @test isa(re1, RanefInfo)
    @test re1.cnames == ["(Intercept)", "days"]
    @test first(re1.levels) == "S308"
    @test first(re1.ranef) ≈ -0.081830 atol = 1e-5
    @test first(re1.stddev) ≈ 0.140354 atol = 1e-5
    f = Figure()
    @test f.content == Any[]
    caterpillar!(f, re1)
    @test length(f.content) == 2
    @test isa(first(f.content), Axis)
    @test size(f.layout) == (1, 2)
    tbl = ranefinfotable(re1)
    @test keys(tbl) == (:name, :level, :cmode, :cstddev)
    @test length(tbl.cmode) == length(re1.cnames) * length(re1.levels)
end
