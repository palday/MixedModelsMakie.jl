f = facetregression(data, :y, :x, :g)
@test save(joinpath(OUTDIR, "facetregression_default.png"), f)

f = facetregression(data, :y, :x, :g; orderby=:slope)
@test save(joinpath(OUTDIR, "facetregression_orderby_slope.png"), f)

f = facetregression(data, :y, :x, :g; orderby=:intercept)
@test save(joinpath(OUTDIR, "facetregression_orderby_intercept.png"), f)

f = facetregression(data, :y, :x, :g; orderby=:slope, rev=true)
@test save(joinpath(OUTDIR, "facetregression_orderby_slope_rev.png"), f)

f = facetregression(data, :y, :x, :g; bank45=false)
@test save(joinpath(OUTDIR, "facetregression_nobank.png"), f)

f = facetregression(data, :y, :x, :g; layout=(5, 4))
@test save(joinpath(OUTDIR, "facetregression_manual_layout.png"), f)

# Partial layout: fix one dimension, auto-compute the other
f = facetregression(data, :y, :x, :g; layout=(2, nothing))
@test save(joinpath(OUTDIR, "facetregression_layout_rows_only.png"), f)

f = facetregression(data, :y, :x, :g; layout=(nothing, 5))
@test save(joinpath(OUTDIR, "facetregression_layout_cols_only.png"), f)

# Undersized explicit layout: never enlarged, warns about the dropped groups
@test_logs((:warn, r"only fits 4 of 20 groups; 16 trailing"),
           begin
               f = facetregression(data, :y, :x, :g; layout=(2, 2))
               @test save(joinpath(OUTDIR, "facetregression_layout_undersized.png"), f)
           end)

# string column names accepted
f = facetregression(data, "y", "x", "g")
@test save(joinpath(OUTDIR, "facetregression_stringnames.png"), f)

# Mutating form into a GridPosition
let fig = Figure()
    facetregression!(fig[1, 1], data, :y, :x, :g)
    @test save(joinpath(OUTDIR, "facetregression_gridpos.png"), fig)
end

@testset "facetregressiontable" begin
    ft = facetregressiontable(data, :y, :x, :g)
    @test names(ft) == ["group", "n", "intercept", "slope"]
    @test nrow(ft) == length(unique(data.g))
    @test all(==(5), ft.n)

    # :none preserves first-encountered order, not alphabetical
    @test ft.group == unique(data.g)

    # :slope sorts increasing
    ft_slope = facetregressiontable(data, :y, :x, :g; orderby=:slope)
    @test issorted(ft_slope.slope)

    # :intercept sorts increasing
    ft_intercept = facetregressiontable(data, :y, :x, :g; orderby=:intercept)
    @test issorted(ft_intercept.intercept)

    # rev reverses whatever ordering orderby produced, including :none
    @test facetregressiontable(data, :y, :x, :g; rev=true).group == reverse(ft.group)
    ft_slope_rev = facetregressiontable(data, :y, :x, :g; orderby=:slope, rev=true)
    @test issorted(ft_slope_rev.slope; rev=true)

    # matches a direct simplelinreg computation for one group
    lab = first(data.g)
    idx = findall(==(lab), data.g)
    a, b = simplelinreg(data.x[idx], data.y[idx])
    row = only(filter(:group => ==(lab), ft))
    @test row.intercept ≈ a
    @test row.slope ≈ b
end

@test_throws ArgumentError facetregression(data, :y, :x, :g; orderby=:bogus)
@test_throws ArgumentError facetregressiontable(data, :y, :x, :g; orderby=:bogus)

# --- LinearMixedModel-based methods ---

f = facetregression(m1, :days)
@test save(joinpath(OUTDIR, "facetregression_model_default.png"), f)

f = facetregression(m1, :days, :subj)
@test save(joinpath(OUTDIR, "facetregression_model_explicit_group.png"), f)

# m0 has a random intercept only (no random slope on days)
f = facetregression(m0, :days)
@test save(joinpath(OUTDIR, "facetregression_model_no_ranef_slope.png"), f)

f = facetregression(m1, :days; show_fixef=false, show_shrunken=false)
@test save(joinpath(OUTDIR, "facetregression_model_toggles_off.png"), f)

# predictor omitted: m1 has exactly one non-intercept fixed effect (days)
f = facetregression(m1)
@test save(joinpath(OUTDIR, "facetregression_model_autopredictor.png"), f)

f = facetregression(m1; group=:subj)
@test save(joinpath(OUTDIR, "facetregression_model_autopredictor_group.png"), f)

# Mutating forms
let fig = Figure()
    facetregression!(fig[1, 1], m1, :days)
    @test save(joinpath(OUTDIR, "facetregression_model_gridpos.png"), fig)
end
let fig = Figure()
    facetregression!(fig[1, 1], m1)
    @test save(joinpath(OUTDIR, "facetregression_model_autopredictor_gridpos.png"), fig)
end

@testset "facetregressiontable (model)" begin
    ft = facetregressiontable(m1, :days)
    @test names(ft) ==
          ["group", "n", "intercept", "slope", "fixef_intercept", "fixef_slope",
           "shrunken_intercept", "shrunken_slope"]
    @test nrow(ft) == length(m1.reterms[1].levels)

    # :none preserves the model's own level order, not resorted
    @test ft.group == m1.reterms[1].levels

    # population fit is constant across rows and matches fixef(m1)
    fe = fixef(m1)
    @test all(==(fe[1]), ft.fixef_intercept)
    @test all(==(fe[2]), ft.fixef_slope)

    # shrunken == fixef + conditional modes, matched by level order
    re = ranef(m1)[1]
    @test ft.shrunken_intercept .- fe[1] ≈ vec(re[1, :])
    @test ft.shrunken_slope .- fe[2] ≈ vec(re[2, :])

    # m0 has no random slope on days: shrunken slope == population slope everywhere
    ft0 = facetregressiontable(m0, :days)
    @test all(≈(fixef(m0)[2]), ft0.shrunken_slope)

    # omitting predictor matches passing it explicitly
    @test facetregressiontable(m1) == facetregressiontable(m1, :days)

    # rev/orderby behave the same as the table-based method
    ft_slope = facetregressiontable(m1, :days; orderby=:slope)
    @test issorted(ft_slope.slope)
    ft_slope_rev = facetregressiontable(m1, :days; orderby=:slope, rev=true)
    @test issorted(ft_slope_rev.slope; rev=true)
end

@test_throws ArgumentError facetregression(m1, :days, :nonexistent)
@test_throws ArgumentError facetregression(m1, :nonexistent_predictor)
@test_throws ArgumentError facetregression(m2) # more than one non-intercept fixed effect
@test_throws ArgumentError facetregressiontable(m1, :days, :nonexistent)
@test_throws ArgumentError facetregressiontable(m1, :nonexistent_predictor)
@test_throws ArgumentError facetregressiontable(m2)
