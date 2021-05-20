module MixedModelsMakie
    using LinearAlgebra
    using Makie
    using MixedModels

    export
        CoefByGroup,
        RanefInfo,

        caterpillar,
        caterpillar!,
        clevelandaxes!,
        ranefinfo,
        shrinkage,
        shrinkageplot,
        shrinkageplot!,
        simplelinreg

    include("shrinkage.jl")
    include("caterpillar.jl")
    include("xyplot.jl")

end # module
