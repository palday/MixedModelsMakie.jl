@test_logs (:warn, "qqline=:R is a deprecated value, use qqline=:fitrobust instead.") match_mode = :any qqnorm(m1;
                                                                                                               qqline=:R)
f = qqnorm(m1; qqline=:fitrobust)
@test save(joinpath(OUTDIR, "qqnorm_sleepstudy_fitrobust.png"), f)

f = qqplot(Normal(0, m1.σ), m1)
@test save(joinpath(OUTDIR, "qqplot_sleepstud.png"), f)
