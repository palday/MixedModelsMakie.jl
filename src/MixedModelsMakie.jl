module MixedModelsMakie

using BSplineKit
using LinearAlgebra
using DataFrames
using Distributions
using KernelDensity
using Makie
using MixedModels
using PrecompileTools
using Printf
using SpecialFunctions
using StatsBase

export RanefInfo,
       caterpillar,
       caterpillar!,
       clevelandaxes!,
       coefplot,
       coefplot!,
       profiledensity,
       profiledensity!,
       qqcaterpillar,
       qqcaterpillar!,
       ranefinfo,
       ranefinfotable,
       regressionplot,
       regressionplot!,
       ridge2d,
       ridge2d!,
       ridgeplot,
       ridgeplot!,
       shrinkageplot,
       shrinkageplot!,
       simplelinreg,
       splom!,
       splomaxes!,
       zetaplot,
       zetaplot!

# from https://github.com/MakieOrg/Makie.jl/issues/2992
const Indexable = Union{Makie.Figure,Makie.GridLayout,Makie.GridPosition,
                        Makie.GridSubposition}

include("utilities.jl")
include("shrinkage.jl")
include("caterpillar.jl")
include("coefplot.jl")
include("profile.jl")
include("ridge.jl")
include("ridge2d.jl")
include("xyplot.jl")
include("recipes.jl")

@setup_workload begin
    model = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days | subj)),
                MixedModels.dataset(:sleepstudy); progress=false)
    @compile_workload begin
        caterpillar(model)
        coefplot(model)
        qqcaterpillar(model)
        qqnorm(model)
        shrinkageplot(model; ellipse=true)
        # not covered:
        # profile plots
        # bootstrap plots (e.g. ridgeplot)
    end
end

end # module
