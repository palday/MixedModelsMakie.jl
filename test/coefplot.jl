f = coefplot(m1)
@test save(joinpath(OUTDIR, "coef_sleepstudy.png"), f)

for legend_pos in (true, false, :top, :bottom, :left, :right, :axis)
    local f = coefplot(m0, m1zc, m1; show_legend=legend_pos,
                       labels=["intercept", "zerocorr", "full"])
    @test save(joinpath(OUTDIR, "coef_sleepstudy_legend_$(legend_pos).png"), f)
end

@test_throws ArgumentError("Inputs differ in coefficient names") coefplot(m1, mr)

f = coefplot(b1)
@test save(joinpath(OUTDIR, "coef_sleepstudy_boot.png"), f)

f = coefplot(mr)
@test save(joinpath(OUTDIR, "coef_rankdeficient.png"), f)

f = coefplot(br)
@test save(joinpath(OUTDIR, "coef_rankdeficient_boot.png"), f)
