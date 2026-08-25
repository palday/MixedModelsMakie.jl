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

## Ridge 2D Plots

```@docs
ridge2d
```

```@example Coefplot
ridge2d(boot)
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

caterpillar!(Figure(; size=(800,600)), subjre)
```

```@example Caterpillar
caterpillar!(Figure(; size=(800,600)), subjre; orderby=2)
```

```@example Caterpillar
caterpillar!(Figure(; size=(800,600)), subjre; orderby=nothing)
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
qqcaterpillar!(Figure(; size=(400,300)), subjre; cols=[1])
```

```@example Caterpillar
qqcaterpillar!(Figure(; size=(400,300)), subjre; cols=[:days])
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
# The SVG sometimes renders incorrectly
CairoMakie.activate!(; type="png")
sleepstudy = MixedModels.dataset(:sleepstudy)
verbagg = MixedModels.dataset(:verbagg)

fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
shrinkageplot(fm1)
```

```@example Shrinkage
shrinkageplot(fm1; ellipse=true)
```

```@example Shrinkage
shrinkageplot!(Figure(; size=(400,400)), fm1)
```

```@example Shrinkage
gm1 = fit(MixedModel,
          @formula(r2 ~ 1 + anger + gender + btype + situ + (1|subj) + (1+gender|item)),
          verbagg,
          Bernoulli();
          progress=false)
shrinkageplot(gm1, :item)
```

## Diagnostics

We have also provided a few useful plot recipes for common plot types applied to mixed models.
These are especially useful for diagnostics and model checking.

### QQ Plots

The methods for `qqnorm` and `qqplot` are implemented using [Makie recipes](https://makie.juliaplots.org/v0.15.0/recipes.html).
In other words, these are convenience wrappers for calling the relevant plotting methods on `residuals(model)`.

Specify the type of line on the QQ plots with the `qqline` keyword-argument. The default for `qqnorm` is `:fitrobust`, which delivers an R-style line connecting the first and third quartiles. The default for `qqplot` is `:identity`, which plots the line with slope = 1 and intercept = 0. The final possibility is `:fit`, which plots the line of best fit (i.e. regressing the quantiles of the residuals onto the quantiles of the reference distribution).

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

## Profile Plots

*Requires MixedModels 4.14 or above.*

```@example Profile
using CairoMakie
using MixedModels
using MixedModelsMakie

sleepstudy = MixedModels.dataset(:sleepstudy)
fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy; progress=false)
pr1 = profile(fm1)
zetaplot(pr1)
```

```@example Profile
# show zeta on the absolute value scale with coverage intervals
zetaplot(pr1; absv=true)
```

```@example Profile
# show zeta for the variance components
zetaplot(pr1; absv=true, ptyp='θ')
```

```@example Profile
profiledensity(pr1)
```

```@example Profile
profiledensity(pr1; share_y_scale=false)
```

```@example Profile
profiledensity(pr1; ptyp='σ')
```

## UpSet Plots

UpSet plots[^upset] visualize the intersection structure of categorical conditions,
showing which combinations of condition levels co-occur in the data and how many
observations (or grouping-factor levels) fall in each combination.

[^upset]: Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister, H. (2014). UpSet: Visualization of Intersecting Sets. IEEE Transactions on Visualization and Computer Graphics, 20(12), 1983–1992. https://doi.org/10.1109/TVCG.2014.2346248

**Sets** are the individual levels of each categorical predictor (e.g., `"gender: M"`,
`"gender: F"`, `"btype: curse"`).
**Columns** of the combination matrix are the full factorial cells — every combination
of one level per predictor — plus *marginal cells* where all levels of a single
predictor are simultaneously active (showing whether that predictor is within-subjects
or between-subjects). A filled circle means that condition level is active in that
column; a connecting line spans the active conditions within each column. The top bars
show counts per column; the left bars show counts per individual condition level.

A column in which only within-subjects predictors are collapsed will have non-zero
counts; a column where a *between*-subjects predictor is collapsed will be empty
because no single unit can appear in all levels of a between-subjects factor.

```@docs
upsetplot
```

```@docs
upsetplot!
```

```@example UpSet
using CairoMakie
CairoMakie.activate!(; type="svg")
using DataFrames: Not
using MixedModels
using MixedModelsMakie

verbagg = MixedModels.dataset(:verbagg)

gm1 = fit(MixedModel,
          @formula(r2 ~ 1 + anger + gender + btype + situ + (1|subj) + (1+gender|item)),
          verbagg, Bernoulli(); progress=false)

# gender is between-subjects; btype and situ are within-subjects.
# Marginal cells that collapse gender are empty (no subject has both genders),
# while marginals that collapse btype or situ are non-empty.
upsetplot(gm1, :subj, show_empty=false)
```

```@example UpSet
# Observation counts instead of subject counts
upsetplot(gm1, nothing, show_empty=false)
```

```@example UpSet
# Table-based: no model required — non-numeric columns are detected automatically.
# Numeric columns (r2, anger) and explicit exclusions (subj, item) are dropped.
upsetplot(verbagg; cols=Not([:subj, :item]), gf=:subj, show_empty=false)
```

```@example UpSet
# Sort by degree (full factorial cells first, then marginals) rather than by count
upsetplot(gm1, :subj, sortby=:degree, show_empty=false)
```

The same incidence data drawn by `upsetplot` can be pulled out as a table with
[`upsettable`](@ref), for cases where you want to filter, sort, or re-analyze it
directly rather than (or in addition to) plotting it:

```@docs
upsettable
```

```@example UpSet
upsettable(gm1, :subj)
```

## Nesting/Crossing Plots

Grouping factors in a mixed model (e.g. `subj`, `item`, `school`, `class`) can be
**nested** — every level of one occurs with exactly one level of another, as with
students nested within classrooms — or **crossed** — levels of both factors
co-occur relatively freely, as with subjects and items in a repeated-measures
design. `nestingplot` shows this structure for every pair of grouping factors in a
model, using only their level assignments (`ReMat.refs`) — no original data frame
is needed, matching the model-based approach used by [`upsetplot`](@ref).

The plot is laid out like a correlation matrix:
- **diagonal**: each grouping factor's name and number of levels.
- **lower triangle**: a heatmap of the co-occurrence contingency table between the
  row and column factor. Levels are reordered (for display only) to group each
  level with its most common partner, so nested structure appears as a
  block-diagonal pattern and crossed structure appears as a dense or scattered
  rectangle.
- **upper triangle**: a text badge classifying the same pair — one factor nested
  in the other (`A ⊂ B`), `identical` (the two names partition observations the
  same way, e.g. two labels for one grouping factor), or crossed, further split
  into `complete` (every combination of levels is observed) and `partial` (only
  some combinations are, with the observed density reported). This split matters
  in practice: a completely-crossed subject × item design supports estimating
  both by-subject and by-item slopes for every condition, while a partially
  crossed (e.g. Latin-square) design may not.

```@docs
nestingplot
```

```@docs
nestingplot!
```

```@example Nesting
using CairoMakie
CairoMakie.activate!(; type="svg")
using DataFrames
using MixedModels
using MixedModelsMakie

kb07 = MixedModels.dataset(:kb07)
gm2 = fit(MixedModel,
          @formula(rt_trunc ~ 1 + spkr * prec * load +
                              (1 + spkr + prec + load | subj) +
                              (1 + spkr | item)),
          kb07; progress=false)

# subj and item are nearly (but not perfectly) crossed in this design
nestingplot(gm2)
```

As with `upsetplot`, the data underlying `nestingplot` is available as tables.
[`nestingtable`](@ref) gives the raw co-occurrence counts in long format, with one
row per pair of levels — including zero-count rows for combinations that never
co-occur, since those absences are exactly what reveal nesting or partial
crossing:

```@docs
nestingtable
```

```@example Nesting
nestingtable(gm2)
```

Filtering to `count == 0` finds specific missing combinations, e.g. the subjects
who never saw a particular item:

```@example Nesting
filter(:count => iszero, nestingtable(gm2))
```

[`nestingstructure`](@ref) gives the pairwise nested/crossed classification and
density shown in `nestingplot`'s upper-triangle badges — one row per pair of
grouping factors, rather than one row per pair of levels:

```@docs
nestingstructure
```

```@example Nesting
nestingstructure(gm2)
```

## General plots

We also provide a `splom` or scatter-plot matrix plot for data frames with numeric columns (i.e. a matrix of all pairwise plots).
These plots can be used to visualize the joint distribution of, say, the parameter estimates from a simulation.

```@docs
splom!
```

```@example Splom
using CairoMakie
using DataFrames
using LinearAlgebra
using MixedModelsMakie

data = rmul!(randn(100, 3), LowerTriangular([+1 +0 +0;
                                             +1 +1 +0;
                                             -1 -1 +1]))
df = DataFrame(data, [:x, :y, :z])

splom!(Figure(; size=(800, 800)), df)
```

Meanwhile, `splomaxes!` provides a lower-level backend for `splom!`

```@docs
splomaxes!
```

```@example Splom
using Statistics

mat = Array(df)
function pfunc(ax, i, j)
    # note that this references mat from the outer scope!
    # [j, i] because:
    # - i is the row in the figure, which corresponds to the y var in each panel
    # - j is the col in the figure, which corresponds to the x var in each panel
    v = view(mat, :, [j, i])
    scatter!(ax, v; color=(:blue, 0.2))
    cc = cor(eachcol(v)...)
    cc = round(cc; digits=2)
    text!(ax, "r=$(cc)")
    return ax
end
splomaxes!(Figure(; size=(800, 800)),
           names(df), pfunc)
```
