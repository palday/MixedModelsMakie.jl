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
CairoMakie.activate!(type = "svg")
using MixedModels
using MixedModelsMakie
sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy)
subjre = ranefinfo(fm1)[:subj]

caterpillar!(Figure(; resolution=(800,600)), subjre)
caterpillar!(Figure(; resolution=(800,600)), subjre; orderby=2)
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

```@example
using CairoMakie
using MixedModels
using MixedModelsMakie
sleepstudy = MixedModels.dataset(:sleepstudy)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy)
shrink = shrinkage(fm1)[:subj]

shrinkageplot!(Figure(; resolution=(800,600)), shrink)
```
