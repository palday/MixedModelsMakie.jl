module MixedModelsMakie
    using LinearAlgebra
    using Makie
    using MixedModels

    export
        RanefInfo,

        caterpillar,
        caterpillar!,
        clevelandaxes!,
        ranefinfo,
        shrinkageplot,
        simplelinreg

    include("shrinkage.jl")
    include("caterpillar.jl")
    include("xyplot.jl")

end # module
