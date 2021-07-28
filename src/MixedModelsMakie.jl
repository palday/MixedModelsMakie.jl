module MixedModelsMakie
    using LinearAlgebra
    using DataFrames
    using Distributions
    using Makie
    using MixedModels
    using SpecialFunctions
    using Statistics

    export
        RanefInfo,

        caterpillar,
        caterpillar!,
        clevelandaxes!,
        coefplot,
        coefplot!,
        qqcaterpillar,
        qqcaterpillar!,
        ranefinfo,
        shrinkageplot,
        shrinkageplot!,
        simplelinreg

    include("utilities.jl")
    include("shrinkage.jl")
    include("caterpillar.jl")
    include("coefplot.jl")
    include("xyplot.jl")
    include("recipes.jl")

end # module
