"""
    coefplot!(f::Figure, x::MixedModel; conf_level=0.95)
    coefplot!(f::Figure, x::MixedModelBootstrap; conf_level=0.95)

Add a "coefplot" of the fixed-effects and associated confidence intervals to the figure.
"""
function coefplot!(f::Figure, x::Union{MixedModel, MixedModelBootstrap}; conf_level=0.95)
	ci = citable(x, conf_level)
	y = nrow(ci):-1:1
	xvals = ci.estimate
	ax = Axis(f[1, 1])
    scatter!(ax, xvals, y)
    errorbars!(ax, xvals, y,  xvals .- ci.lower, ci.upper .- xvals, direction=:x)
    ax.xlabel = "Estimate and $(conf_level * 100)% confidence interval"
    ax.yticks = (y, ci.coefname)
	ylims!(0, nrow(ci) + 1)
    f
end

_npreds(x::MixedModelBootstrap) = length(first(x.fits).β)
_npreds(x::MixedModel) = length(coefnames(x))

"""
    coefplot(x::MixedModel; conf_level=0.95)
    coefplot(x::MixedModelBootstrap; conf_level=0.95)

Return a `Figure` of a "coefplot" of the fixed-effects and associated confidence intervals.
"""
function coefplot(x::Union{MixedModel, MixedModelBootstrap}; conf_level=0.95)
	coefplot!(Figure(resolution=(640, 75 * _npreds(x))), x)
end

"""
    citable(x::MixedModel, level=0.95)
    citable(x::MixedModelBootstrap, level=0.95)

Return a DataFrame of coefficient names, point estimates and confidence intervals.

`level` specifies the confidence level, e.g. `0.95` for 95% confidence intervals.

For `MixedModels`, the intervals are computed using the standard errors and the asymptotic
Wald approximation (e.g. est±1.96*se for 95% intervals). For `MixedModelBootstrap`, the intervals
are computed using `shortestcovint`, and the point estimate re-computed as the mean of the bootstrap values.

The returned table has the following columns:
- `coefname`: the names of the coefficients
- `estimate`: the point estimates
- `lower`: the lower edge of the confidence interval
- `upper`: the upper edge of the confidence interval

"""
function citable(x::MixedModel, level=0.95)
	# taking from the lower tail
	semultiple = zquantile((1 - level) / 2)
	se = stderror(x)

	return DataFrame(coefname=coefnames(x), estimate=coef(x),
	                 lower=coef(x) + semultiple * se,
	                 upper=coef(x) - semultiple * se)
end

function citable(x::MixedModelBootstrap, level=0.95)
	df = transform!(select!(DataFrame(x.β), Not(:iter)),
		            :coefname => ByRow(string) => :coefname)
	combine(groupby(df, :coefname),
		    :β => mean => :estimate,
		    :β => NamedTuple{(:lower, :upper)} ∘ shortestcovint => [:lower, :upper])

end
