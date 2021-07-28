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

```@docs
coefplot!
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

coefplot(boot; conf_level=0.99)
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

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
subjre = ranefinfo(fm1)[:subj]

caterpillar!(Figure(; resolution=(800,600)), subjre)
```

```@example Caterpillar
caterpillar!(Figure(; resolution=(800,600)), subjre; orderby=2)
```

```@example Caterpillar
caterpillar!(Figure(; resolution=(800,600)), subjre; orderby=nothing)
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

## Diagnostics

We have also provided a few useful plot recipes for common plot types applied to mixed models.
These are especially useful for diagnostics and model checking.
### QQ Plots

```@example Residuals
using CairoMakie
using MixedModels
using MixedModelsMakie

sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
qqnorm(fm1)
```

```@example Residuals
# the residuals should have mean 0
# and standard deviation equal to the residual standard deviation
qqplot(Normal(0, fm1.σ), fm1)
```
