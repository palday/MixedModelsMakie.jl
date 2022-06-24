using CairoMakie
using MixedModels
using MixedModelsMakie
using Random # we don't depend on exact PRNG vals, so no need for StableRNGs
using Test

using MixedModelsMakie: confint_table

const OUTDIR = joinpath(pkgdir(MixedModelsMakie), "test", "output")

@testset "utilities, types and tables" begin
    include("utils_and_types.jl")
end

m1 = fit(MixedModel,
         @formula(1000 / reaction ~ 1 + days + (1 + days | subj)),
         MixedModels.dataset(:sleepstudy))

m2 = fit(MixedModel,
         @formula(rt_trunc ~ 1 + spkr * prec * load +
                             (1 + spkr + prec + load | subj) +
                             (1 + spkr | item)),
         MixedModels.dataset(:kb07))

b1 = parametricbootstrap(MersenneTwister(42), 100, m1)

@testset "[qq]caterpillar" begin
    f = caterpillar(m1)
    save(joinpath(OUTDIR, "cat_sleepstudy.png"), f)

    f = caterpillar(m2, :subj)
    save(joinpath(OUTDIR, "cat_kb07_subj.png"), f)

    f = caterpillar(m2, :item)
    save(joinpath(OUTDIR, "cat_kb07_item.png"), f)

    f = qqcaterpillar(m1)
    save(joinpath(OUTDIR, "qqcat_sleepstudy.png"), f)

    f = qqcaterpillar(m2, :subj)
    save(joinpath(OUTDIR, "qqcat_kb07_subj.png"), f)

    f = qqcaterpillar(m2, :item)
    save(joinpath(OUTDIR, "qqcat_kb07_item.png"), f)
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
end
