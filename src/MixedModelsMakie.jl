module MixedModelsMakie
    using AbstractPlotting
    using LinearAlgebra
    using MixedModels

    export
        CoefByGroup,
        RanefInfo,

        caterpillar,
        caterpillar!,
        ranefinfo,
        shrinkage,
        shrinkageplot

    include("shrinkage.jl")
    include("caterpillar.jl")

end
