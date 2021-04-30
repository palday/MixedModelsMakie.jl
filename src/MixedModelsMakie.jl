module MixedModelsMakie
    using AbstractPlotting
    using MixedModels

    export
        CoefByGroup,

        shrinkage,
        shrinkageplot

    include("shrinkage.jl")

end
