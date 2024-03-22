"""
    coefplot(x::Union{MixedModel,MixedModelBootstrap}; kwargs...)::Figure
    coefplot!(fig::$(Indexable), x::Union{MixedModel,MixedModelBootstrap}; 
              kwargs...)
    coefplot!(ax::Axis, Union{MixedModel,MixedModelBootstrap};
              conf_level=0.95, vline_at_zero=true, show_intercept=true, attributes...)

Create a coefficient plot of the fixed-effects and associated confidence intervals.

`attributes` are passed onto `scatter!` and `errorbars!`.

The mutating methods return the original object.
"""
function coefplot(x::Union{MixedModel,MixedModelBootstrap}; show_intercept=true, kwargs...)
    # need to guarantee a min height of 150
    fig = Figure(; size=(640, max(150, 75 * _npreds(x; show_intercept))))
    coefplot!(fig, x; kwargs...)
    return fig
end

"""$(@doc coefplot)"""
function coefplot!(fig::Indexable, x::Union{MixedModel,MixedModelBootstrap}; kwargs...)
    ax = Axis(fig[1, 1])
    coefplot!(ax, x; kwargs...)
    return fig
end

"""$(@doc coefplot)"""
function coefplot!(ax::Axis, x::Union{MixedModel,MixedModelBootstrap};
                   conf_level=0.95,
                   vline_at_zero=true,
                   show_intercept=true,
                   attributes...)
    ci = confint_table(x, conf_level)
    show_intercept || filter!(:coefname => !=("(Intercept)"), ci)
    y = nrow(ci):-1:1
    xvals = ci.estimate
    xlabel = @sprintf "Estimate and %g%% confidence interval" conf_level * 100

    attributes = merge((; xlabel), attributes)
    ax.xlabel = attributes.xlabel

    scatter!(ax, xvals, y; attributes...)
    errorbars!(ax, xvals, y, xvals .- ci.lower, ci.upper .- xvals;
               direction=:x, attributes...)
    vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)

    reset_limits!(ax)
    cn = _coefnames(x; show_intercept)
    nticks = _npreds(x; show_intercept)
    ax.yticks = (nticks:-1:1, cn)
    ylims!(ax, 0, nticks + 1)
    return ax
end
