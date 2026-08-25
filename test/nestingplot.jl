using MixedModelsMakie: _cooccurrence, _nesting_relationship, _nesting_data

@testset "_nesting_relationship" begin
    # perfect nesting: each column (B) has exactly one nonzero row (A)
    nested = [1 1 0 0
              0 0 1 0
              0 0 0 1]
    rel = _nesting_relationship(nested)
    @test rel.rel === :b_nested_in_a
    @test rel.density ≈ 4 / 12

    # transpose: A nested in B
    rel_t = _nesting_relationship(permutedims(nested))
    @test rel_t.rel === :a_nested_in_b

    # identical partition (permutation table): both directions hold
    ident = [3 0 0
             0 5 0
             0 0 2]
    rel_i = _nesting_relationship(ident)
    @test rel_i.rel === :identical

    # complete crossing: every cell nonzero
    complete = [1 1 1
                1 1 1]
    rel_c = _nesting_relationship(complete)
    @test rel_c.rel === :complete_crossing
    @test rel_c.density == 1.0

    # partial crossing: dense but not every combination observed
    partial = [1 1 0
               1 1 1]
    rel_p = _nesting_relationship(partial)
    @test rel_p.rel === :partial_crossing
    @test rel_p.density ≈ 5 / 6
end

@testset "_cooccurrence" begin
    tab = _cooccurrence([1, 1, 2, 2, 2], 2, [1, 2, 1, 2, 2], 2)
    @test tab == [1 1
                  1 2]
end

# m2: kb07, subj + item are (nearly completely) crossed grouping factors
f = nestingplot(m2)
@test save(joinpath(OUTDIR, "nesting_kb07.png"), f)

info = MixedModelsMakie._nesting_data(m2)
@test info.names == ("subj", "item")
@test info.relationships[(2, 1)].rel === :partial_crossing

# g1: verbagg (GLMM), subj + item grouping
f = nestingplot(g1)
@test save(joinpath(OUTDIR, "nesting_verbagg.png"), f)

# gfs restricts/orders the grouping factors shown
f = nestingplot(m2; gfs=[:item, :subj])
@test save(joinpath(OUTDIR, "nesting_kb07_gfs.png"), f)

# Mutating form into a GridPosition
let f = Figure()
    nestingplot!(f[1, 1], m2)
    @test save(joinpath(OUTDIR, "nesting_kb07_gridpos.png"), f)
end

# Synthetic fully-nested design: class nested in school, unique class ids
let rng = MersenneTwister(1)
    n_school, n_class, n_student = 4, 3, 10
    school = String[]
    class = String[]
    y = Float64[]
    for s in 1:n_school, c in 1:n_class, _ in 1:n_student
        push!(school, "school$s")
        push!(class, "school$(s)_class$(c)")
        push!(y, randn(rng))
    end
    data = (; y, school, class)
    mnest = fit(MixedModel, @formula(y ~ 1 + (1 | school) + (1 | class)), data;
                progress)
    info = MixedModelsMakie._nesting_data(mnest)
    rel = only(values(info.relationships))
    @test rel.rel === :b_nested_in_a || rel.rel === :a_nested_in_b

    f = nestingplot(mnest)
    @test save(joinpath(OUTDIR, "nesting_synthetic_nested.png"), f)
end

# Errors
@test_throws ArgumentError nestingplot(m1) # only one distinct grouping factor
@test_throws ArgumentError nestingplot(m2; gfs=[:subj, :nonexistent])
