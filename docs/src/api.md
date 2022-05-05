```@meta
CurrentModule = MixedModelsMakie
DocTestSetup = quote
    using MixedModelsMakie
end
DocTestFilters = [r"([a-z]*) => \1", r"getfield\(.*##[0-9]+#[0-9]+"]
```

# MixedModelsMakie.jl API

## Coefficient Plots

```@docs
coefplot
```

```@example Coefplot
using CairoMakie
using MixedModels
using MixedModelsMakie
using Random

verbagg = MixedModels.dataset(:verbagg)

gm1 = fit(MixedModel,
          @formula(r2 ~ 1 + anger + gender + btype + situ + (1|subj) + (1|item)),
          verbagg,
          Bernoulli();
          progress=false)

coefplot(gm1)
```

```@example Coefplot
sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
boot = parametricbootstrap(MersenneTwister(42), 1000, fm1)

coefplot(boot; conf_level=0.999, title="Custom Title")
```

## Ridge Plots

```@docs
ridgeplot
```

```@example Coefplot
ridgeplot(boot)
```


## Random effects and group-level predictions

### Caterpillar Plots

```@docs
RanefInfo
```

```@docs
ranefinfo
```

```@docs
caterpillar
```

```@docs
caterpillar!
```

```@example Caterpillar
using CairoMakie
CairoMakie.activate!(type = "svg")
using MixedModels
using MixedModelsMakie
sleepstudy = MixedModels.dataset(:sleepstudy)
verbagg = MixedModels.dataset(:verbagg)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
gm0 = fit(MixedModel,
          @formula(r2 ~ 1 + anger + gender + btype + situ + (1|subj) + (1|item)),
          verbagg,
          Bernoulli();
          progress=false)

subjre = ranefinfo(fm1)[:subj]

caterpillar!(Figure(; resolution=(800,600)), subjre)
```

```@example Caterpillar
caterpillar!(Figure(; resolution=(800,600)), subjre; orderby=2)
```

```@example Caterpillar
caterpillar!(Figure(; resolution=(800,600)), subjre; orderby=nothing)
```

```@example Caterpillar
caterpillar(gm0, :item)
```

```@docs
qqcaterpillar
```

```@docs
qqcaterpillar!
```

```@example Caterpillar
qqcaterpillar(fm1)
```

```@example Caterpillar
qqcaterpillar(gm0, :item)
```

```@example Caterpillar
qqcaterpillar!(Figure(; resolution=(400,300)), subjre; cols=[1])
```

```@example Caterpillar
qqcaterpillar!(Figure(; resolution=(400,300)), subjre; cols=[:days])
```

### Shrinkage Plots

```@docs
shrinkageplot
shrinkageplot!
```

```@example Shrinkage
using CairoMakie
using MixedModels
using MixedModelsMakie
sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
shrinkageplot(fm1)
```

```@example Shrinkage
shrinkageplot!(Figure(; resolution=(400,400)), fm1)
```

```@example Shrinkage
gm1 = fit(MixedModel,
          @formula(r2 ~ 1 + anger + gender + btype + situ + (1|subj) + (1+gender|item)),
          verbagg,
          Bernoulli();
          progress=false)
shrinkageplot(gm1)
```

## Diagnostics

We have also provided a few useful plot recipes for common plot types applied to mixed models.
These are especially useful for diagnostics and model checking.

### QQ Plots

The methods for `qqnorm` and `qqplot` are implemented using [Makie recipes](https://makie.juliaplots.org/v0.15.0/recipes.html).
In other words, these are convenience wrappers for calling the relevant plotting methods on `residuals(model)`.

Specify the type of line on the QQ plots with the `qqline` keyword-argument. The default for `qqnorm` is `:fitrobust`, which delivers an R-style line connecting the first and third quartiles. The default for `qqplot` is `:identity`, which plots the line with slope = 1 and intercept = 0. The final possiblity is `:fit`, which plots the line of best fit (i.e. regressing the quantiles of the residuals onto the quantiles of the reference distribution).

The reference distribution for `qqnorm` is the standard normal, which differs from [the behavior in previous versions of Makie](https://github.com/JuliaPlots/Makie.jl/pull/1277).

!!! compat
    The [options and associated names for the `qqline` keyword argument](https://makie.juliaplots.org/v0.16/examples/plotting_functions/qqplot/index.html) changed in [Makie 0.16.3](https://github.com/JuliaPlots/Makie.jl/pull/1563) (and were broken in Makie 0.16.0-0.16.2). The equivalent to `qqline=:R` is `qqline=:fitrobust`. `qqline=:R` will be supported for backwards compatibility only until the next breaking release.
```@example Residuals
using CairoMakie
using MixedModels
using MixedModelsMakie

sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
qqnorm(fm1; qqline=:fitrobust)
```

```@example Residuals
# the residuals should have mean 0
# and standard deviation equal to the residual standard deviation
qqplot(Normal(0, fm1.σ), fm1)
```

## General plots

We also provide a `splom` or scatter-plot matrix plot for data frames with numeric columns (i.e. a matrix of all pairwise plots).
These plots can be used to visualize the joint distribution of, say, the parameter estimates from a simulation.

```@example Splom
using CairoMakie
using DataFrames
using LinearAlgebra
using MixedModelsMakie

splom!(
    Figure(; resolution=(800, 800)),
    DataFrame(rmul!(randn(100, 3), LowerTriangular([1 0 0;1 1 0;-1 -1 1])), [:x, :y, :z]),
)
```
