module MixedModelsMakie
    using LinearAlgebra
    using Makie
    using MixedModels
    using SpecialFunctions

    export
        RanefInfo,

        caterpillar,
        caterpillar!,
        clevelandaxes!,
        qqcaterpillar,
        qqcaterpillar!,
        ranefinfo,
        shrinkageplot,
        shrinkageplot!,
        simplelinreg

    include("utilities.jl")
    include("shrinkage.jl")
    include("caterpillar.jl")
    include("xyplot.jl")

end # module
