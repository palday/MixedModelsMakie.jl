"""
    coefplot(x::Union{MixedModel,MixedModelBootstrap}; kwargs...)::Figure
    coefplot!(fig::$(Indexable), x::Union{MixedModel,MixedModelBootstrap};
              kwargs...)
    coefplot!(ax::Axis, Union{MixedModel,MixedModelBootstrap};
              conf_level=0.95, vline_at_zero=true, show_intercept=true,
              scatter_attributes=(;),
              errorbars_attributes=(;),
              attributes...)

Create a coefficient plot of the fixed-effects and associated confidence intervals.

Inestimable coefficients (coefficients removed by pivoting in the rank deficient case) are excluded.

`attributes` are passed onto both `scatter!` and `errorbars!`, while
`scatter_attributes` and `errorbars_attributes` are passed only onto `scatter!` and
`errorbars!`, respectively. (Starting with Makie 0.21, unsupported attributes for a
given plottype are no longer silently ignored, so it's necessary so separate out the
attributes that are only valid for a single plot type.)

The mutating methods return the original object.

!!! note
    Inestimable coefficients (coefficients removed by pivoting in the rank deficient case)
    are excluded.
"""
function coefplot(x::Union{MixedModel,MixedModelBootstrap}; show_intercept=true, kwargs...)
    # need to guarantee a min height of 150
    fig = Figure(; size=(640, max(150, 75 * _npreds(x; show_intercept))))
    coefplot!(fig, x; show_intercept, kwargs...)
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
                   scatter_attributes=(;),
                   errorbars_attributes=(;),
                   attributes...)
    ci = confint_table(x, conf_level; show_intercept)
    y = nrow(ci):-1:1
    xvals = ci.estimate
    xlabel = @sprintf "Estimate and %g%% confidence interval" conf_level * 100

    ax.xlabel = xlabel

    scatter!(ax, xvals, y; attributes..., scatter_attributes...)
    errorbars!(ax, xvals, y, xvals .- ci.lower, ci.upper .- xvals;
               direction=:x, attributes..., errorbars_attributes...)
    vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)

    reset_limits!(ax)
    cn = _coefnames(x; show_intercept)
    nticks = _npreds(x; show_intercept)
    ax.yticks = (nticks:-1:1, cn)
    ylims!(ax, 0, nticks + 1)
    return ax
end
