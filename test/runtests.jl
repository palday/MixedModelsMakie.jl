include("setup_tests.jl")

@testset ExtendedTestSet "MixedModelsMakie.jl" begin
    @testset "Aqua" begin
        # we can't check for unbound type parameters
        # because we actually need one at one point for _same_family()
        Aqua.test_all(MixedModelsMakie;
                      ambiguities=false,
                      unbound_args=false,
                      piracies=(; treat_as_own=[MixedModel]))
    end

    @testset "utilities, types and tables" include("utils_and_types.jl")

    @testset "[qq]caterpillar" include("caterpillar.jl")

    @testset "clevelandaxes" include("clevelandaxes.jl")

    @testset "facetregression" include("facetregression.jl")

    @testset "coefplot" include("coefplot.jl")

    @testset "recipes" include("recipes.jl")

    @testset "ridgeplot" include("ridgeplot.jl")

    @testset "ridge2d" include("ridge2d.jl")

    @testset "shrinkageplot" include("shrinkageplot.jl")

    @testset "splom!" include("splom.jl")

    @testset "profile" include("profile.jl")

    @testset "upsetplot" include("upsetplot.jl")

    @testset "nestingplot" include("nestingplot.jl")
end
