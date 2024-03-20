@test_throws(ArgumentError("No parameters x found."),
             ridge2d(b1; ptype=:x))
@test_throws(ArgumentError("Only 1 ρ-parameter found: 2D plots require at least 2."),
             ridge2d(b1; ptype=:ρ))
@test save(joinpath(OUTDIR, "ridge2d_beta.png"), ridge2d(b1))
@test save(joinpath(OUTDIR, "ridge2d_sigma.png"), ridge2d(b1; ptype=:σ))
