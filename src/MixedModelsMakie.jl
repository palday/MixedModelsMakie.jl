module MixedModelsMakie
using AlgebraOfGraphics
using BSplineKit
using LinearAlgebra
using DataFrames
using Distributions
using KernelDensity
using Makie
using MixedModels
using Printf
using SpecialFunctions
using StatsBase

export RanefInfo,
       caterpillar,
       caterpillar!,
       clevelandaxes!,
       coefplot,
       coefplot!,
       deltadensity,
       qqcaterpillar,
       qqcaterpillar!,
       ranefinfo,
       ranefinfotable,
       ridgeplot,
       ridgeplot!,
       shrinkageplot,
       shrinkageplot!,
       simplelinreg,
       splom!,
       zetaplot,
       zetatraceplot

include("utilities.jl")
include("shrinkage.jl")
include("caterpillar.jl")
include("coefplot.jl")
include("profile.jl")
include("ridge.jl")
include("xyplot.jl")
include("recipes.jl")

end # module
