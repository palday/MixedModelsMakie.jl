f = coefplot(m1)
@test save(joinpath(OUTDIR, "coef_sleepstudy.png"), f)

f = coefplot(b1)
@test save(joinpath(OUTDIR, "coef_sleepstudy_boot.png"), f)
