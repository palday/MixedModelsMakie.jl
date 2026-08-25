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

@testset "nestingtable" begin
    nt = nestingtable(m2)
    @test names(nt) == ["factor_a", "level_a", "factor_b", "level_b", "count"]
    # one row per (level of subj) x (level of item) combination
    @test nrow(nt) == length(m2.reterms[1].levels) * length(m2.reterms[2].levels)
    # co-occurrence counts sum to the number of observations
    @test sum(nt.count) == nobs(m2)

    # gfs restricts which factors are compared (factor_a/factor_b assignment
    # is positional, not gfs order, so only check the *set* of factors shown)
    nt2 = nestingtable(m2; gfs=[:item, :subj])
    @test Set(vcat(nt2.factor_a, nt2.factor_b)) == Set(["subj", "item"])

    # errors mirror nestingplot
    @test_throws ArgumentError nestingtable(m1)
    @test_throws ArgumentError nestingtable(m2; gfs=[:subj, :nonexistent])

    # exact nesting is visible as zero-count rows for non-matching pairs
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
        ntnest = nestingtable(mnest)
        @test nrow(ntnest) == n_school * (n_school * n_class)
        @test count(iszero, ntnest.count) > 0
        @test count(!iszero, ntnest.count) == n_school * n_class
    end
end

@testset "nestingstructure" begin
    ns = nestingstructure(m2)
    @test names(ns) ==
          ["factor_a", "factor_b", "n_levels_a", "n_levels_b", "relationship",
           "density"]
    # one row per pair of grouping factors
    @test nrow(ns) == 1
    @test only(ns.relationship) === :partial_crossing
    @test only(ns.density) ≈ 1789 / 1792

    # gfs restricts which factors are compared
    ns2 = nestingstructure(m2; gfs=[:item, :subj])
    @test Set(vcat(ns2.factor_a, ns2.factor_b)) == Set(["subj", "item"])

    # errors mirror nestingplot
    @test_throws ArgumentError nestingstructure(m1)
    @test_throws ArgumentError nestingstructure(m2; gfs=[:subj, :nonexistent])

    # exact nesting is correctly classified
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
        nsnest = nestingstructure(mnest)
        rel = only(nsnest.relationship)
        @test rel === :b_nested_in_a || rel === :a_nested_in_b
    end
end
