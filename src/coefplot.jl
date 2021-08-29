_npreds(x::MixedModelBootstrap) = length(first(x.fits).β)
_npreds(x::MixedModel) = length(coefnames(x))

_coefnames(x) = coefnames(x)
_coefnames(x::MixedModelBootstrap) = collect(string.(propertynames(first(x.fits).β)))

"""
    coefplot(x::MixedModel; conf_level=0.95, vline_at_zero=true, attributes...)
    coefplot(x::MixedModelBootstrap; conf_level=0.95, vline_at_zero=true, attributes...)

Create a coefficient plot of the fixed-effects and associated confidence intervals.

!!! note
    This functionality is implemented using [Makie recipes](https://makie.juliaplots.org/v0.15.0/recipes.html)
    and thus there are also additional auto-generated methods for `coefplot` and `coefplot!` that may be useful
    when constructing more complex figures.
"""
function coefplot(
    x::Union{MixedModel,MixedModelBootstrap};
    conf_level=0.95,
    vline_at_zero=true,
    attributes...,
)
    fig = Figure(; resolution=(640, 75 * _npreds(x)))
    ax = Axis(fig[1, 1])
    pl = coefplot!(ax, x; conf_level, vline_at_zero, attributes...)
    return Makie.FigureAxisPlot(fig, ax, pl)
end

@recipe(CoefPlot, x) do scene
    return Attributes(; conf_level=0.95, vline_at_zero=true)
end

function Makie.plot!(ax::Axis, P::Type{<:CoefPlot}, allattrs::Makie.Attributes, x)
    plot = Makie.plot!(ax.scene, P, allattrs, x)

    if haskey(allattrs, :title)
        ax.title = allattrs.title[]
    end
    if haskey(allattrs, :xlabel)
        ax.xlabel = allattrs.xlabel[]
    else
        ax.xlabel = @sprintf "Estimate and %g%% confidence interval" (
            allattrs.conf_level[] * 100
        )
    end
    if haskey(allattrs, :ylabel)
        ax.ylabel = allattrs.ylabel[]
    end
    reset_limits!(ax)
    cn = _coefnames(x)
    nticks = _npreds(x)
    ax.yticks = (nticks:-1:1, cn)
    ylims!(ax, 0, nticks + 1)
    allattrs.vline_at_zero[] && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)
    return plot
end

function Makie.plot!(plot::CoefPlot{<:Tuple{Union{MixedModel,MixedModelBootstrap}}})
    model_or_boot, conf_level = plot[1][], plot.conf_level[]
    ci = confint_table(model_or_boot, conf_level)
    y = nrow(ci):-1:1
    xvals = ci.estimate
    scatter!(plot, xvals, y)
    errorbars!(plot, xvals, y, xvals .- ci.lower, ci.upper .- xvals; direction=:x)

    return plot
end

# TODO: it would be useful to have a ridge plot of by-coefficient densities for a richer plot of the bootstrap

"""
    confint_table(x::StatsBase.StatisticalModel, level=0.95)
    confint_table(x::MixedModelBootstrap, level=0.95)

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
function confint_table(x::StatsBase.StatisticalModel, level=0.95)
    # taking from the lower tail
    semultiple = zquantile((1 - level) / 2)
    se = stderror(x)

    return DataFrame(;
        coefname=coefnames(x),
        estimate=coef(x),
        lower=coef(x) + semultiple * se,
        upper=coef(x) - semultiple * se,
    )
end

function confint_table(x::MixedModelBootstrap, level=0.95)
    df = transform!(
        select!(DataFrame(x.β), Not(:iter)), :coefname => ByRow(string) => :coefname
    )
    return combine(
        groupby(df, :coefname),
        :β => mean => :estimate,
        :β => NamedTuple{(:lower, :upper)} ∘ shortestcovint => [:lower, :upper],
    )
end
