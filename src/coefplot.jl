"""
    coefplot(xs::Union{MixedModel,MixedModelBootstrap}...; kwargs...)::Figure
    coefplot!(fig::$(Indexable), xs::Union{MixedModel,MixedModelBootstrap}...;
              show_legend=length(xs) > 1, legend_attributes=(;), kwargs...)
    coefplot!(ax::Axis, xs::Union{MixedModel,MixedModelBootstrap}...;
              conf_level=0.95, vline_at_zero=true, show_intercept=true,
              scatter_attributes=(;),
              errorbars_attributes=(;),
              show_legend=length(xs) > 1,
              legend_attributes=(;),
              labels=string.(1:length(xs)),
              attributes...)

Create a coefficient plot of the fixed-effects and associated confidence intervals.
When multiple models are supplied, they are overlaid on the same axes for comparison;
all models must share the same coefficient names.

`attributes` are passed onto both `scatter!` and `errorbars!`, while
`scatter_attributes` and `errorbars_attributes` are passed only onto `scatter!` and
`errorbars!`, respectively. (Starting with Makie 0.21, unsupported attributes for a
given plottype are no longer silently ignored, so it's necessary to separate out the
attributes that are only valid for a single plot type.)

`labels` controls the legend entry for each model (defaults to `"1"`, `"2"`, ...).

`show_legend` controls placement of the legend. Accepted values:
- `false`: no legend
- `true` or `:bottom`: horizontal legend below the axis (default for multi-model figures)
- `:top`: horizontal legend above the axis
- `:left`: vertical legend to the left of the axis
- `:right`: vertical legend to the right of the axis
- `:axis`: legend embedded inside the axis via `axislegend`

`legend_attributes` is a named tuple of keyword arguments forwarded to the Makie
`Legend` (or `axislegend`) constructor.

The mutating methods return the original object.

!!! note
    Inestimable coefficients (coefficients removed by pivoting in the rank deficient case)
    are excluded. When multiple models are provided, they must have the same
    coefficient names _after_ dropping inestimable coefficients.
"""
function coefplot(xs::Union{MixedModel,MixedModelBootstrap}...;
                  show_intercept=true,
                  show_legend=length(xs) > 1,
                  kwargs...)
    width = 640
    # need to guarantee a min height of 150
    height = max(150, 75 * _npreds(first(xs); show_intercept))

    width += 50 * (show_legend in (:left, :right))
    height += 50 * (show_legend in (true, :top, :bottom))

    fig = Figure(; size=(width, height))
    coefplot!(fig, xs...; show_intercept, show_legend, kwargs...)
    return fig
end

"""$(@doc coefplot)"""
function coefplot!(fig::Indexable, xs::Union{MixedModel,MixedModelBootstrap}...;
                   show_legend=true, legend_attributes=(;), kwargs...)
    ax = Axis(fig[1, 1])
    kwargs = _extract_title!(ax, kwargs)
    axis_legend = show_legend === :axis
    legend_attributes = merge((; merge=true, unique=true), legend_attributes)
    if axis_legend
        show_legend = false
    end
    coefplot!(ax, xs...; legend_attributes, show_legend=axis_legend, kwargs...)

    _place_legend!(fig, ax, show_legend; legend_attributes...)
    return fig
end

"""$(@doc coefplot)"""
function coefplot!(ax::Axis, xs::Union{MixedModel,MixedModelBootstrap}...;
                   conf_level=0.95,
                   vline_at_zero=true,
                   show_intercept=true,
                   scatter_attributes=(;),
                   errorbars_attributes=(;),
                   show_legend=length(xs) > 1,
                   legend_attributes=(;),
                   labels=string.(1:length(xs)),
                   attributes...)
    x = first(xs)
    cn = _coefnames(x; show_intercept)
    nticks = _npreds(x; show_intercept)
    all(_coefnames(m; show_intercept) == cn for m in xs) ||
        throw(ArgumentError("Inputs differ in coefficient names"))

    if length(xs) == 1
        color = get(attributes, :color, :black)
        attributes = merge((; color=color), attributes)
    end

    for (x, label) in zip(xs, labels)
        ci = confint_table(x, conf_level; show_intercept)
        y = nrow(ci):-1:1
        xvals = ci.estimate
        scatter!(ax, xvals, y; attributes..., scatter_attributes..., label)
        errorbars!(ax, xvals, y, xvals .- ci.lower, ci.upper .- xvals;
                   direction=:x, attributes..., errorbars_attributes...,
                   label)
    end

    vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)

    # this is the axis method, so we only have the axislegend
    if show_legend
        legend_attributes = merge((; merge=true, unique=true), legend_attributes)
        axislegend(ax; legend_attributes...)
    end

    reset_limits!(ax)
    xlabel = @sprintf "Estimate and %g%% confidence interval" conf_level * 100
    ax.xlabel = xlabel
    ax.yticks = (nticks:-1:1, cn)
    ylims!(ax, 0, nticks + 1)
    return ax
end

# function g(; kwargs...)
#     if :a in keys(kwargs)
#         a = kwargs[:a]
#         kwargs = Base.pairs(NamedTuple((k => v for (k, v) in kwargs if k != :a)))
#     end

#     @info "" kwargs
# end
