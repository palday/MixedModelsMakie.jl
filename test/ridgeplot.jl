f = ridgeplot(b1)
@test save(joinpath(OUTDIR, "ridge_sleepstudy.png"), f)

f = ridgeplot(br; color=(:blue, 0.3), errorbars_attributes=(; whiskerwidth=15))
@test save(joinpath(OUTDIR, "ridge_rankdeficient.png"), f)

for legend_pos in (true, false, :top, :bottom, :left, :right, :axis)
    local f = ridgeplot(b0, b1; show_legend=legend_pos, labels=["intercept-only", "full"])
    @test save(joinpath(OUTDIR, "ridge_sleepstudy_legend_$(legend_pos).png"), f)
end

@test_throws ArgumentError("Inputs differ in coefficient names") ridgeplot(b1, br)
