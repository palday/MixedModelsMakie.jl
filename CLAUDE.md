# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Run the full test suite:
```julia
julia --project -e'using Pkg; Pkg.test()'
```

Run a single test file (e.g., caterpillar tests):
```julia
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/setup_tests.jl"); include("test/test_caterpillar.jl")'
```

Format code (uses YAS style via JuliaFormatter):
```julia
julia --project -e 'using JuliaFormatter; format(".")'
```

Build documentation:
```julia
julia --project=docs docs/make.jl
```

## Architecture

MixedModelsMakie.jl provides Makie-based visualizations for mixed-effects models from MixedModels.jl. The package is WIP; breaking changes are expected before 1.0, and minor version bumps can be breaking.

### Three-tier plotting API

Every plot type follows this pattern:
1. `function(model; kwargs...)` — creates a `Figure`, delegates to mutating form
2. `function!(fig::Indexable, model; kwargs...)` — creates `Axis`, delegates to axis form
3. `function!(ax::Axis, model; kwargs...)` — does the actual drawing

`Indexable` is a union of `Figure`, `GridLayout`, and `GridPosition` for flexible embedding.

### Plot types by file

- [src/caterpillar.jl](src/caterpillar.jl) — horizontal error bar plots of random effects (conditional means ± CI), plus QQ-scale variant. Uses `RanefInfo` struct to cache extracted ranef data.
- [src/coefplot.jl](src/coefplot.jl) — fixed-effect coefficient plots with CI; accepts both `MixedModel` and `MixedModelBootstrap`.
- [src/shrinkage.jl](src/shrinkage.jl) — scatter-plot matrix of unshrunken vs. shrunken RE estimates with correlation ellipses. `_ranef()` temporarily mutates model θ to compute unshrunken values — not thread-safe.
- [src/ridge.jl](src/ridge.jl) — ridge density plots of bootstrap parameter distributions; optionally overlays coefplot.
- [src/ridge2d.jl](src/ridge2d.jl) — bivariate scatter + contour plots for bootstrap parameter pairs.
- [src/profile.jl](src/profile.jl) — profile likelihood ζ plots with BSplineKit interpolation; density plots of profile ζ values.
- [src/xyplot.jl](src/xyplot.jl) — general-purpose utilities: `clevelandaxes!` (multi-panel grid layout), `splom!`/`splomaxes!` (scatter-plot matrix accepting custom panel functions), `simplelinreg`.
- [src/utilities.jl](src/utilities.jl) — shared helpers: `confint_table` (extracts CIs to DataFrame), `_coefnames`, `ppoints`, `zquantile`, `_extract_title!`.
- [src/upset.jl](src/upset.jl) — UpSet plots of categorical predictor intersection structure. Model-based API recovers predictor levels from the contrast matrix (`_obs_predictor_cols`); table-based API auto-detects non-numeric columns. `gf` is a positional argument (like `shrinkageplot`).
- [src/nesting.jl](src/nesting.jl) — nesting/crossing matrix for a model's grouping factors, using only `ReMat.refs`/`.levels` (no data needed). Correlation-matrix-style layout: diagonal names each factor, lower triangle heatmaps the co-occurrence table (reordered to reveal block-diagonal nesting), upper triangle badges the pairwise classification (nested/identical/complete or partial crossing).
- [src/recipes.jl](src/recipes.jl) — `convert_arguments` extensions for Makie's `QQNorm`/`QQPlot` to accept `MixedModel` objects.

### Testing

Tests live in [test/](test/); each plot type has its own file (e.g., `test_caterpillar.jl`). [test/setup_tests.jl](test/setup_tests.jl) fits shared models (`m1`, `m1zc`, `m2`, `g1`) and bootstrap samples (`b1`, `br`) used across all test files. Reference images for visual regression are stored in [test/output/](test/output/) and compared via Percy.io in CI.

### Code style

YAS style (enforced by `style.yml` CI workflow via JuliaFormatter + reviewdog).
