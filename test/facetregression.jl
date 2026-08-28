f = facetregression(data, :y, :x, :g)
@test save(joinpath(OUTDIR, "facetregression_default.png"), f)

f = facetregression(data, :y, :x, :g; orderby=:slope)
@test save(joinpath(OUTDIR, "facetregression_orderby_slope.png"), f)

f = facetregression(data, :y, :x, :g; orderby=:intercept)
@test save(joinpath(OUTDIR, "facetregression_orderby_intercept.png"), f)

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
