module MixedModelsMakie
    using LinearAlgebra
    using Makie
    using MixedModels
    using SpecialFunctions

    export
        RanefInfo,

        CaterpillarPlot,
        caterpillarplot,
        caterpillarplot!,
        clevelandaxes!,
        #QQCaterpillarPlot,
        #qqcaterpillar,
        #qqcaterpillar!,
        ranefinfo,
        ShrinkagePlot,
        shrinkageplot,
        shrinkageplot!,
        simplelinreg

    include("utilities.jl")
    include("shrinkage.jl")
    include("caterpillar.jl")
    include("xyplot.jl")

end # module
