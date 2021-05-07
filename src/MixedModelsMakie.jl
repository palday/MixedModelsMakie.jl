module MixedModelsMakie
    using AbstractPlotting
    using LinearAlgebra
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
        simplelinreg

    include("shrinkage.jl")
    include("caterpillar.jl")
    include("xyplot.jl")

end
