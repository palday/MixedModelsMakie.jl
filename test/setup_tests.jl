using Aqua
using CairoMakie
using DataFrames
using MixedModels
using MixedModelsMakie
using Random # we don't depend on exact PRNG vals, so no need for StableRNGs
using Statistics
using Suppressor
using Test
using TestSetExtensions

using MixedModelsMakie: confint_table

const OUTDIR = joinpath(pkgdir(MixedModelsMakie), "test", "output")
const progress = false

function save(path, obj, args...; kwargs...)
    isfile(path) && rm(path)
    Makie.save(path, obj, args...; kwargs...)
    return isfile(path)
end

m1 = fit(MixedModel,
         @formula(1000 / reaction ~ 1 + days + (1 + days | subj)),
         MixedModels.dataset(:sleepstudy); progress)

m1zc = fit(MixedModel,
           @formula(1000 / reaction ~ 1 + days + zerocorr(1 + days | subj)),
           MixedModels.dataset(:sleepstudy); progress)

m2 = fit(MixedModel,
         @formula(rt_trunc ~ 1 + spkr * prec * load +
                             (1 + spkr + prec + load | subj) +
                             (1 + spkr | item)),
         MixedModels.dataset(:kb07); progress)

b1 = parametricbootstrap(MersenneTwister(42), 500, m1; progress,
                         optsum_overrides=(; ftol_rel=1e-6))
g1 = fit(MixedModel,
         @formula(r2 ~ 1 + anger + gender + btype + situ +
                       (1 | subj) + (1 + gender | item)),
         MixedModels.dataset(:verbagg),
         Bernoulli(); progress)

rng = MersenneTwister(0)
x = rand(rng, 100)
a = rand(rng, 100)
data = (; x, x2=1.5 .* x, a, y=rand(rng, [0, 1], 100), g=repeat('A':'T', 5))
mr = @suppress fit(MixedModel, @formula(y ~ a + x + x2 + (1 | g)), data; progress)
br = @suppress parametricbootstrap(MersenneTwister(42), 500, mr; progress,
                                   optsum_overrides=(; ftol_rel=1e-6))
