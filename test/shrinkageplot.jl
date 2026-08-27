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

@test_throws(ArgumentError("You only specified a single column. You need at least two."),
             shrinkageplot(m2; ellipse=true, cols=["spkr: old"]))

@test_throws(ArgumentError("Grouping variable (\"subj\") only has a single predictor associated with it (\"(Intercept)\"). You need at least two."),
             shrinkageplot(m2int; ellipse=false))

@test_throws(ArgumentError("Grouping variable (\"item\") only has a single predictor associated with it (\"(Intercept)\"). You need at least two."),
             shrinkageplot(m2int, :item))

f = shrinkageplot(m2; ellipse=true, ellipse_scale=2)
@test save(joinpath(OUTDIR, "shrinkage_kb07_subj_ellipse_scaled.png"), f)

f = shrinkageplot(g1, :item)
@test save(joinpath(OUTDIR, "shrinkage_verbagg.png"), f)

f = shrinkageplot(g1, :item; ellipse=true, n_ellipse=2)
@test save(joinpath(OUTDIR, "shrinkage_verbagg_ellipse.png"), f)

f = shrinkageplot(m1; labels=true)
@test save(joinpath(OUTDIR, "shrinkage_sleepstudy_labels.png"), f)

f = shrinkageplot(m2, :subj; labels=["S030", "S031"])
@test save(joinpath(OUTDIR, "shrinkage_kb07_subj_labels_subset.png"), f)

@test_throws(ArgumentError("Specified columns not found in random effects: [\"nonexistent\"]"),
             shrinkageplot(m1; labels=["nonexistent"]))
