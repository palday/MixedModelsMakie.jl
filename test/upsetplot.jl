kb07data = MixedModels.dataset(:kb07)

# m2: kb07, spkr*prec*load categorical fixed effects, subj + item grouping
f = upsetplot(m2)
@test save(joinpath(OUTDIR, "upset_kb07_subj.png"), f)

f = upsetplot(m2, :item)
@test save(joinpath(OUTDIR, "upset_kb07_item.png"), f)

# g1: verbagg (GLMM), gender + btype + situ categorical fixed effects
f = upsetplot(g1)
@test save(joinpath(OUTDIR, "upset_verbagg_subj.png"), f)

# Mutating form into a GridPosition
let f = Figure()
    upsetplot!(f[1, 1], m2)
    @test save(joinpath(OUTDIR, "upset_kb07_gridpos.png"), f)
end

# Observation counts (gf=nothing)
f = upsetplot(m2, nothing)
@test save(joinpath(OUTDIR, "upset_kb07_obs.png"), f)

# Table-based: auto-exclude numeric rt_trunc; restrict to condition columns
f = upsetplot(kb07data; cols=Not([:subj, :item]))
@test save(joinpath(OUTDIR, "upset_kb07_table_obs.png"), f)

# Table-based with gf — count unique subjects per cell
f = upsetplot(kb07data; cols=Not([:subj, :item]), gf=:subj)
@test save(joinpath(OUTDIR, "upset_kb07_table_gf.png"), f)

# Table-based mutating form
let fig = Figure()
    upsetplot!(fig[1, 1], kb07data; cols=Not([:subj, :item]))
    @test save(joinpath(OUTDIR, "upset_kb07_table_gridpos.png"), fig)
end

# Error: no categorical columns left after filtering (only numerics)
@test_throws ArgumentError upsetplot(DataFrame(; x=rand(10), y=rand(10)))

# Error: unknown grouping factor
@test_throws ArgumentError upsetplot(m2, :nonexistent)

# Error: model has no categorical fixed effects (m1 uses only continuous `days`)
@test_throws ArgumentError upsetplot(m1)

# Error: too many factorial cells (high-cardinality columns included)
@test_throws ArgumentError upsetplot(kb07data; cols=Not(:rt_trunc))

# Warning: high-cardinality columns detected
@test_logs((:warn, r"subj has 56 distinct levels"),
           (:warn, r"item has 32 distinct levels"),
           match_mode = :any,
           @test_throws ArgumentError upsetplot(kb07data; cols=Not(:rt_trunc)))
