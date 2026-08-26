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

# --- Layout options ---

# Vertical orientation
f = upsetplot(m2; orientation=:vertical)
@test save(joinpath(OUTDIR, "upset_kb07_vertical.png"), f)

# Vertical with intersection on the right
f = upsetplot(m2; orientation=:vertical, intersection_pos=:right)
@test save(joinpath(OUTDIR, "upset_kb07_vertical_intright.png"), f)

# Vertical with set-size bars on the bottom (x-axis labels move to top)
f = upsetplot(m2; orientation=:vertical, setsize_pos=:bottom)
@test save(joinpath(OUTDIR, "upset_kb07_vertical_setbot.png"), f)

# Horizontal with intersection on the bottom
f = upsetplot(m2; intersection_pos=:bottom)
@test save(joinpath(OUTDIR, "upset_kb07_intbottom.png"), f)

# Horizontal with set-size bars on the right
f = upsetplot(m2; setsize_pos=:right)
@test save(joinpath(OUTDIR, "upset_kb07_setright.png"), f)

# Horizontal with both flipped
f = upsetplot(m2; intersection_pos=:bottom, setsize_pos=:right)
@test save(joinpath(OUTDIR, "upset_kb07_intbottom_setright.png"), f)

# No set-size bar chart
f = upsetplot(m2; show_setsize=false)
@test save(joinpath(OUTDIR, "upset_kb07_nosetsize.png"), f)

# Vertical without set-size
f = upsetplot(m2; orientation=:vertical, show_setsize=false)
@test save(joinpath(OUTDIR, "upset_kb07_vertical_nosetsize.png"), f)

# Error: invalid orientation
@test_throws ArgumentError upsetplot(m2; orientation=:diagonal)

# Error: wrong position for orientation
@test_throws ArgumentError upsetplot(m2; intersection_pos=:left)
@test_throws ArgumentError upsetplot(m2; orientation=:vertical, intersection_pos=:top)
@test_throws ArgumentError upsetplot(m2; setsize_pos=:top)
@test_throws ArgumentError upsetplot(m2; orientation=:vertical, setsize_pos=:left)

@testset "upsettable" begin
    ut = upsettable(g1)
    @test names(ut)[1:3] == ["cell", "degree", "count"]
    @test Symbol("gender: M") in propertynames(ut)
    @test all(x -> x isa Bool, ut[!, Symbol("gender: M")])
    # counts and cell order match the data underlying upsetplot itself
    @test ut.count == MixedModelsMakie._upset_data(g1).cell_counts

    # gf=nothing counts observations instead of grouping-factor levels
    ut_obs = upsettable(g1, nothing)
    @test ut_obs.count != ut.count

    # table-based method matches the model-based one in shape
    ut_tab = upsettable(kb07data; cols=Not([:subj, :item]))
    @test names(ut_tab)[1:3] == ["cell", "degree", "count"]

    @test_throws ArgumentError upsettable(m1) # no categorical fixed effects
end
