using CairoMakie, MixedModels, MixedModelsMakie, Random
progress = true

m1 = fit(MixedModel,
         @formula(1000 / reaction ~ 1 + days + (1 + days | subj)),
         MixedModels.dataset(:sleepstudy); progress)


b1 = parametricbootstrap(MersenneTwister(42), 500, m1; progress,
                         optsum_overrides=(;ftol_rel=1e-6))

@test_throws(ArgumentError("No parameters x found."),
             ridge2d(b1; ptype=:x))
@test_throws(ArgumentError("Only 1 ρ-paramater found: 2D plots require at least 2."),
             ridge2d(b1; ptype=:ρ))
ridge2d(b1)
ridge2d(b1; ptype=:σ)
