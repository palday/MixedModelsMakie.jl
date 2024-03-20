f = caterpillar(m1; vline_at_zero=true)
@test save(joinpath(OUTDIR, "cat_sleepstudy.png"), f)

f = caterpillar(m2, :subj)
@test save(joinpath(OUTDIR, "cat_kb07_subj.png"), f)

f = caterpillar(m2, :item)
@test save(joinpath(OUTDIR, "cat_kb07_item.png"), f)

f = caterpillar(m2, :subj; cols=[:("load: yes"), :("prec: maintain")], orderby=2)
@test save(joinpath(OUTDIR, "cat_kb07_subj_ordered_cols.png"), f)

@test_throws ArgumentError caterpillar(m2, :subj; cols=[:("load: no")])

f = caterpillar(g1)
@test save(joinpath(OUTDIR, "cat_verbagg.png"), f)

f = qqcaterpillar(m1; vline_at_zero=true)
@test save(joinpath(OUTDIR, "qqcat_sleepstudy.png"), f)

f = qqcaterpillar(m2, :subj)
@test save(joinpath(OUTDIR, "qqcat_kb07_subj.png"), f)

f = qqcaterpillar(m2, :item)
@test save(joinpath(OUTDIR, "qqcat_kb07_item.png"), f)

f = qqcaterpillar(g1)
@test save(joinpath(OUTDIR, "qqcat_verbagg.png"), f)

let f = Figure(; size=(1000, 600))
    gl = f[1, 1] = GridLayout()
    re = ranefinfo(m2)
    qqcaterpillar!(gl, re[:item])
    Label(gl[end + 1, :], "Item"; font=:bold)
    gl = f[1, 2] = GridLayout()
    qqcaterpillar!(gl, re[:subj])
    Label(gl[end + 1, :], "Subject"; font=:bold)
    Label(f[0, :], "Conditional Modes")
    colsize!(f.layout, 1, Auto(0.5))
    save(joinpath(OUTDIR, "qqcat_kb07_joint.png"), f)
    f
end
