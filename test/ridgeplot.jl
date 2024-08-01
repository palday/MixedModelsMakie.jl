f = ridgeplot(b1)
@test save(joinpath(OUTDIR, "ridge_sleepstudy.png"), f)

f = ridgeplot(br; color=(:blue, 0.3), errorbars_attributes=(; whiskerwidth=15))
@test save(joinpath(OUTDIR, "ridge_rankdeficient.png"), f)
