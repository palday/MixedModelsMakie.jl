```@meta
CurrentModule = MixedModelsMakie
DocTestSetup = quote
    using MixedModelsMakie
end
DocTestFilters = [r"([a-z]*) => \1", r"getfield\(.*##[0-9]+#[0-9]+"]
```

# MixedModelsMakie.jl API

## Caterpillar Plots

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

```@example
using CairoMakie
using MixedModels
using MixedModelsMakie
sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy)

caterpillar(fm1)
```

## Shrinkage Plots

```@docs
CoefByGroup
```

```@docs
shrinkage
```

```@docs
shrinkageplot
```

```@docs
shrinkageplot!
```

```@docs
shrinkage2d!
```

```@example
using CairoMakie
using MixedModels
using MixedModelsMakie
sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy)

shrinkageplot(fm1)
```
