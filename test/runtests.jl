using CairoMakie
using DataFrames
using MixedModels
using MixedModelsMakie
using Random # we don't depend on exact PRNG vals, so no need for StableRNGs
using Statistics
using Test

using MixedModelsMakie: confint_table

const OUTDIR = joinpath(pkgdir(MixedModelsMakie), "test", "output")
const progress = false

@testset "utilities, types and tables" begin
    include("utils_and_types.jl")
end

m1 = fit(MixedModel,
         @formula(1000 / reaction ~ 1 + days + (1 + days | subj)),
         MixedModels.dataset(:sleepstudy); progress)

m2 = fit(MixedModel,
         @formula(rt_trunc ~ 1 + spkr * prec * load +
                             (1 + spkr + prec + load | subj) +
                             (1 + spkr | item)),
         MixedModels.dataset(:kb07); progress)

b1 = parametricbootstrap(MersenneTwister(42), 100, m1)

g1 = fit(MixedModel,
         @formula(r2 ~ 1 + anger + gender + btype + situ +
                       (1 | subj) + (1 + gender | item)),
         MixedModels.dataset(:verbagg),
         Bernoulli(); progress)

@testset "[qq]caterpillar" begin
    f = caterpillar(m1)
    save(joinpath(OUTDIR, "cat_sleepstudy.png"), f)

    f = caterpillar(m2, :subj)
    save(joinpath(OUTDIR, "cat_kb07_subj.png"), f)

    f = caterpillar(m2, :item)
    save(joinpath(OUTDIR, "cat_kb07_item.png"), f)

    f = caterpillar(g1)
    save(joinpath(OUTDIR, "cat_verbagg.png"), f)

    f = qqcaterpillar(m1)
    save(joinpath(OUTDIR, "qqcat_sleepstudy.png"), f)

    f = qqcaterpillar(m2, :subj)
    save(joinpath(OUTDIR, "qqcat_kb07_subj.png"), f)

    f = qqcaterpillar(m2, :item)
    save(joinpath(OUTDIR, "qqcat_kb07_item.png"), f)

    f = qqcaterpillar(g1)
    save(joinpath(OUTDIR, "qqcat_verbagg.png"), f)
end

@testset "clevelandaxes" begin
    f = clevelandaxes!(Figure(), ["S$(lpad(i, 2))" for i in 1:16], (4, 4))
    n = 12
    for i in 1:4, j in 1:4
        x = randn(MersenneTwister(i), n)
        y = randn(MersenneTwister(j), n)
        scatter!(f[i, j], x, y)
        text!(f[i, j], 1.9, -1.9;
              text="[$i, $j]", align=(:center, :center), textsize=14)
    end
    save(joinpath(OUTDIR, "clevelandaxes.png"), f)
end

@testset "coefplot" begin
    f = coefplot(m1)
    save(joinpath(OUTDIR, "coef_sleepstudy.png"), f)

    f = coefplot(b1)
    save(joinpath(OUTDIR, "coef_sleepstudy_boot.png"), f)
end

@testset "recipes" begin
    @test_logs (:warn, "qqline=:R is a deprecated value, use qqline=:fitrobust instead.") match_mode = :any qqnorm(m1;
                                                                                                                   qqline=:R)
    f = qqnorm(m1; qqline=:fitrobust)
    save(joinpath(OUTDIR, "qqnorm_sleepstudy_fitrobust.png"), f)

    f = qqplot(Normal(0, m1.σ), m1)
    save(joinpath(OUTDIR, "qqplot_sleepstud.png"), f)
end

@testset "ridgeplot" begin
    f = ridgeplot(b1)
    save(joinpath(OUTDIR, "ridge_sleepstudy.png"), f)
end

@testset "shrinkageplot" begin
    f = shrinkageplot(m1)
    save(joinpath(OUTDIR, "shrinkage_sleepstudy.png"), f)

    f = shrinkageplot(m2, :item)
    save(joinpath(OUTDIR, "shrinkage_kb07_item.png"), f)

    f = shrinkageplot(m2, :subj)
    save(joinpath(OUTDIR, "shrinkage_kb07_subj.png"), f)

    f = shrinkageplot(g1, :item)
    save(joinpath(OUTDIR, "shrinkage_verbagg.png"), f)
end

@testset "splom!" begin
    df = DataFrame(MixedModels.dataset(:mmec))
    splof = @test_logs (:info,
                        r"Ignoring 3 non-numeric columns") splom!(Figure(), df)
    save(joinpath(OUTDIR, "splom_mmec.png"), splof)
end
