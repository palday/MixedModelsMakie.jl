# MixedModelsMakie

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable Docs][docs-stable-img]][docs-stable-url]
[![Dev Docs][docs-dev-img]][docs-dev-url]
[![Codecov](https://codecov.io/gh/palday/MixedModelsMakie.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/palday/MixedModelsMakie.jl)
[![DOI](https://zenodo.org/badge/337082315.svg)](https://zenodo.org/badge/latestdoi/337082315)
[![Code Style: YAS](https://img.shields.io/badge/code%20style-yas-1fdcb2.svg)](https://github.com/jrevels/YASGuide)

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://palday.github.io/MixedModelsMakie.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://palday.github.io/MixedModelsMakie.jl/stable

`MixedModelsMakie.jl` is a Julia package providing plotting capabilities for models fit with [MixedModels.jl](https://juliastats.org/MixedModels.jl/stable/).
Plotting is performed using the [Makie ecoysystem](https://makie.juliaplots.org/stable/), with the interface defined using `Makie.jl` ([previously `AbstractPlotting.jl`](https://discourse.julialang.org/t/ann-makie-v-0-13/61522)) for compatibility across all Makie backends.

Note that the functionality here is currently early development and so breaking changes are expected as we refine the interface.
Following [SemVer](https://semver.org/), these minor releases before 1.0 can introduce these breaking changes.
The release of 1.0 will be indication that we believe the interface is reasonably stable.
Minor refinements of graphical displays without changing the API are considered non-breaking.
