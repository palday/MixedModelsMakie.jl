f = shrinkageplot(m1)
@test save(joinpath(OUTDIR, "shrinkage_sleepstudy.png"), f)

f = shrinkageplot(m2, :item)
@test save(joinpath(OUTDIR, "shrinkage_kb07_item.png"), f)

f = shrinkageplot(m2, :subj)
@test save(joinpath(OUTDIR, "shrinkage_kb07_subj.png"), f)

f = shrinkageplot(m2; ellipse=true)
@test save(joinpath(OUTDIR, "shrinkage_kb07_subj_ellipse.png"), f)

f = shrinkageplot(m2; ellipse=true, cols=["spkr: old", "prec: maintain", "(Intercept)"])
@test save(joinpath(OUTDIR, "shrinkage_kb07_subj_cols.png"), f)

@test_throws(ArgumentError("At least two columns must be specified."),
             shrinkageplot(m2; ellipse=true, cols=["spkr: old"]))

f = shrinkageplot(m2; ellipse=true, ellipse_scale=2)
@test save(joinpath(OUTDIR, "shrinkage_kb07_subj_ellipse_scaled.png"), f)

f = shrinkageplot(g1, :item)
@test save(joinpath(OUTDIR, "shrinkage_verbagg.png"), f)

f = shrinkageplot(g1, :item; ellipse=true, n_ellipse=2)
@test save(joinpath(OUTDIR, "shrinkage_verbagg_ellipse.png"), f)
