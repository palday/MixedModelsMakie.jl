f = ridgeplot(b1)
@test save(joinpath(OUTDIR, "ridge_sleepstudy.png"), f)
