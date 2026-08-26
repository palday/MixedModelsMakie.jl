@testset "utilities" begin
    ppts = MixedModelsMakie.ppoints(64)
    @test length(ppts) == 64
    @test first(ppts) ≈ inv(128)

    @test MixedModelsMakie.zquantile(0.025) ≈ -1.96 atol = 0.005
    @test MixedModelsMakie.zquantile(0.975) ≈ 1.96 atol = 0.005
    @test MixedModelsMakie.zquantile(0.50) ≈ 0
end

@testset "_histcurve" begin
    curve = MixedModelsMakie._histcurve(zeros(50); bins=1)
    @test curve isa MixedModelsMakie._StepCurve
    @test length(curve.x) == length(curve.density)
    # closed at 0 on both ends, and the single bin spans the whole range
    @test first(curve.density) == 0
    @test last(curve.density) == 0
    @test maximum(curve.density) == 50

    curve2 = MixedModelsMakie._histcurve(zeros(50); bins=5)
    @test maximum(curve2.density) == 50
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

@testset "confint_table" begin
    wald = confint_table(m1_speed, 0.68)
    bsamp = parametricbootstrap(MersenneTwister(42), 1000, m1_speed; progress)
    boot = confint_table(bsamp, 0.68)

    @test wald.coefname == boot.coefname
    @test wald.estimate ≈ boot.estimate rtol = 0.05
    @test wald.lower ≈ boot.lower rtol = 0.05
    @test wald.upper ≈ boot.upper rtol = 0.05

    @test all(splat(isapprox),
              zip(MixedModelsMakie.confint_table(mr).estimate, fixef(mr)))

    @test fixefnames(mr) == MixedModelsMakie.confint_table(mr).coefname
    @test fixefnames(mr) == MixedModelsMakie.confint_table(br).coefname

    sigma = confint_table(b1, 0.68; ptype=:σ)
    @test sort(sigma.coefname) ==
          ["residual", "subj: (Intercept)", "subj: days"]

    rho = confint_table(b1, 0.68; ptype=:ρ)
    @test rho.coefname == ["subj: (Intercept), days"]

    theta = confint_table(b1, 0.68; ptype=:θ)
    @test sort(theta.coefname) == ["θ1", "θ2", "θ3"]

    @test_throws ArgumentError confint_table(b1; ptype=:σs)

    sigma_subj = confint_table(b2, 0.68; ptype=:σ, group=:subj)
    @test sort(sigma_subj.coefname) ==
          ["subj: (Intercept)", "subj: load: yes", "subj: prec: maintain",
           "subj: spkr: old"]

    sigma_item = confint_table(b2, 0.68; ptype=:σ, group=:item)
    @test sort(sigma_item.coefname) == ["item: (Intercept)", "item: spkr: old"]

    rho_item = confint_table(b2, 0.68; ptype=:ρ, group=:item)
    @test rho_item.coefname == ["item: (Intercept), spkr: old"]

    @test_throws ArgumentError confint_table(b2; ptype=:β, group=:subj)
    @test_throws ArgumentError confint_table(b2; ptype=:θ, group=:subj)
    @test_throws ArgumentError confint_table(b2; ptype=:σ, group=:nope)

    @test confint_table(b1, 0.68; ptype=:sigma) == sigma
    @test confint_table(b1, 0.68; ptype=:rho) == rho
    @test confint_table(b1, 0.68; ptype=:theta) == theta
end

@testset "ranefinfo" begin
    reinfo = ranefinfo(m1_speed)
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
