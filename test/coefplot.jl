f = coefplot(m1)
@test save(joinpath(OUTDIR, "coef_sleepstudy.png"), f)

f = coefplot(b1)
@test save(joinpath(OUTDIR, "coef_sleepstudy_boot.png"), f)

f = coefplot(mr)
@test save(joinpath(OUTDIR, "coef_rankdeficient.png"), f)

f = coefplot(br)
@test save(joinpath(OUTDIR, "coef_rankdeficient_boot.png"), f)
