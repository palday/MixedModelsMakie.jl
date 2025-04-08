var documenterSearchIndex = {"docs":
[{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"CurrentModule = MixedModelsMakie\nDocTestSetup = quote\n    using MixedModelsMakie\nend\nDocTestFilters = [r\"([a-z]*) => \\1\", r\"getfield\\(.*##[0-9]+#[0-9]+\"]","category":"page"},{"location":"api/#MixedModelsMakie.jl-API","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"","category":"section"},{"location":"api/#Caterpillar-Plots","page":"MixedModelsMakie.jl API","title":"Caterpillar Plots","text":"","category":"section"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"RanefInfo","category":"page"},{"location":"api/#MixedModelsMakie.RanefInfo","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.RanefInfo","text":"RanefInfo\n\nInformation on random effects conditional modes/means, variances, etc.\n\nUsed for creating caterpillar plots.\n\nnote: Note\nThis functionality may be moved upstream into MixedModels.jl in the near future.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"ranefinfo","category":"page"},{"location":"api/#MixedModelsMakie.ranefinfo","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.ranefinfo","text":"ranefinfo(m::LinearMixedModel)\n\nReturn a NamedTuple{fnames(m), NTuple(k, RanefInfo)} from model m\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"caterpillar","category":"page"},{"location":"api/#MixedModelsMakie.caterpillar","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.caterpillar","text":"caterpillar(m::LinearMixedModel, gf::Symbol)\n\nReturns a Figure of a \"caterpillar plot\" of the random-effects means and prediction intervals\n\nA \"caterpillar plot\" is a horizontal error-bar plot of conditional means and standard deviations of the random effects.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"caterpillar!","category":"page"},{"location":"api/#MixedModelsMakie.caterpillar!","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.caterpillar!","text":"caterpillar!(f::Figure, r::RanefInfo; orderby=1)\n\nAdd Axes of a caterpillar plot from r to f.\n\nThe order of the levels on the vertical axes is increasing orderby column of r.ranef, usually the (Intercept) random effects. Setting orderby=nothing will disable sorting, i.e. return the levels in the order they are stored in.\n\nnote: Note\nEven when not sorting the levels, they might have already been sorted during model matrix construction. If you want impose a particular ordering on the levels, then you must sort the relevant fields in the RanefInfo object before calling caterpillar!.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"using CairoMakie\nCairoMakie.activate!(type = \"svg\")\nusing MixedModels\nusing MixedModelsMakie\nsleepstudy = MixedModels.dataset(:sleepstudy)\n\nfm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy)\nsubjre = ranefinfo(fm1)[:subj]\n\ncaterpillar!(Figure(; resolution=(800,600)), subjre)","category":"page"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"caterpillar!(Figure(; resolution=(800,600)), subjre; orderby=2)","category":"page"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"caterpillar!(Figure(; resolution=(800,600)), subjre; orderby=nothing)","category":"page"},{"location":"api/#Shrinkage-Plots","page":"MixedModelsMakie.jl API","title":"Shrinkage Plots","text":"","category":"section"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"shrinkageplot","category":"page"},{"location":"api/#MixedModelsMakie.shrinkageplot","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.shrinkageplot","text":"shrinkageplot(m::LinearMixedModel, gf::Symbol=first(fnames(m)), θref)\n\nReturn a scatter-plot matrix of the conditional means, b, of the random effects for grouping factor gf.\n\nTwo sets of conditional means are plotted: those at the estimated parameter values and those at θref. The default θref results in Λ being a very large multiple of the identity.  The corresponding conditional means can be regarded as unpenalized.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"MixedModelsMakie.jl API","title":"MixedModelsMakie.jl API","text":"using CairoMakie\nusing MixedModels\nusing MixedModelsMakie\nsleepstudy = MixedModels.dataset(:sleepstudy)\n\nfm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy)\nshrinkageplot(fm1)","category":"page"},{"location":"#MixedModelsMakie.jl-Documentation","page":"MixedModelsMakie.jl Documentation","title":"MixedModelsMakie.jl Documentation","text":"","category":"section"},{"location":"","page":"MixedModelsMakie.jl Documentation","title":"MixedModelsMakie.jl Documentation","text":"CurrentModule = MixedModelsMakie","category":"page"},{"location":"","page":"MixedModelsMakie.jl Documentation","title":"MixedModelsMakie.jl Documentation","text":"MixedModelsMakie.jl is a Julia package providing plotting capabilities for models fit with MixedModels.jl.","category":"page"},{"location":"","page":"MixedModelsMakie.jl Documentation","title":"MixedModelsMakie.jl Documentation","text":"Note that the functionality here is currently early alpha development and so breaking changes are expected as we refine the interface. Following SemVer, these minor releases before 1.0 can introduce these breaking changes. The release of 1.0 will be indication that we believe the interface is reasonably stable.","category":"page"},{"location":"","page":"MixedModelsMakie.jl Documentation","title":"MixedModelsMakie.jl Documentation","text":"Pages = [\n        \"api.md\",\n]\nDepth = 1","category":"page"}]
}
