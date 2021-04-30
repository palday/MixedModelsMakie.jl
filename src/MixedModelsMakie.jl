module MixedModelsMakie
    using CairoMakie
    using MixedModels

    export
        CoefByGroup,

        shrinkage,
        shrinkageplot

    include("shrinkage.jl")

end
