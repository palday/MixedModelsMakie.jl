f = ridgeplot(b1)
@test save(joinpath(OUTDIR, "ridge_sleepstudy.png"), f)

f = ridgeplot(br)
@test save(joinpath(OUTDIR, "ridge_rankdeficient.png"), f)
