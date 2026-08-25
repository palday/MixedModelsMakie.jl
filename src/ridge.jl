"""
    ridgeplot(x::Union{MixedModel,MixedModelBootstrap}...; kwargs...)::Figure
    ridgeplot!(fig::$(Indexable), x::MixedModelBootstrap...;
               show_legend=length(xs) > 1, legend_attributes=(;), kwargs...)
    ridgeplot!(ax::Axis, x::MixedModelBootstrap...;
               conf_level=0.95, vline_at_zero=true, show_intercept=true,
               scatter_attributes=(;),
               errorbars_attributes=(;),
               band_attributes=(;),
               lines_attributes=(;),
               show_legend=length(xs) > 1,
               legend_attributes=(;),
               labels=string.(1:length(xs)),
               attributes...)

Create a ridge plot for the bootstrap samples of the fixed effects.
When multiple bootstrap objects are supplied, they are overlaid on the same axes for
comparison; all inputs must share the same coefficient names.

Densities are normalized so that the maximum density is always 1.

The highest density interval corresponding to `conf_level` is marked with a bar at the bottom of each density.
Setting `conf_level=missing` removes the markings for the highest density interval.

`attributes` are passed onto [`coefplot`](@ref), `band!` and `lines!`.
`scatter_attributes` and `errorbars_attributes` are passed only onto [`coefplot`](@ref).
`band_attributes` and `lines_attributes` are passed only onto `band!` and
`lines!`, respectively.
(Starting with Makie 0.21, unsupported attributes for a
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
    are excluded.
"""
function ridgeplot(xs::MixedModelBootstrap...;
                   show_intercept=true,
                   show_legend=length(xs) > 1,
                   kwargs...)
    width = 640
    height = max(200, 100 * _npreds(first(xs); show_intercept))

    width += 50 * (show_legend in (:left, :right))
    height += 50 * (show_legend in (true, :top, :bottom))

    fig = Figure(; size=(width, height))
    return ridgeplot!(fig, xs...; show_intercept, show_legend, kwargs...)
end

"""$(@doc ridgeplot)"""
function ridgeplot!(fig::Indexable, xs::MixedModelBootstrap...; show_legend=true,
                    legend_attributes=(;), kwargs...)
    ax = Axis(fig[1, 1])
    kwargs = _extract_title!(ax, kwargs)
    axis_legend = show_legend === :axis
    legend_attributes = merge((; merge=true, unique=true), legend_attributes)
    if axis_legend
        show_legend = false
    end
    ridgeplot!(ax, xs...; legend_attributes, show_legend=axis_legend, kwargs...)

    _place_legend!(fig, ax, show_legend; legend_attributes...)
    return fig
end

"""
    _color(s::Symbol)
    _color(p::Pair)

Extract the color part out of either a color name or a `(color, alpha)` pair.
"""
_color(s) = s
_color(p::Pair) = first(p)

"""$(@doc ridgeplot)"""
function ridgeplot!(ax::Axis, xs::MixedModelBootstrap...;
                    conf_level=0.95,
                    vline_at_zero=true,
                    show_intercept=true,
                    scatter_attributes=(;),
                    errorbars_attributes=(;),
                    band_attributes=(;),
                    lines_attributes=(;),
                    show_legend=length(xs) > 1,
                    legend_attributes=(;),
                    labels=string.(1:length(xs)),
                    attributes...)
    x = first(xs)
    cn = _coefnames(x; show_intercept)
    all(_coefnames(m; show_intercept) == cn for m in xs) ||
        throw(ArgumentError("Inputs differ in coefficient names"))

    xlabel = if !ismissing(conf_level)
        @sprintf "Normalized bootstrap density and %g%% confidence interval" (conf_level *
                                                                              100)
    else
        "Normalized bootstrap density"
    end

    if length(xs) == 1
        color = get(attributes, :color, :black)
        attributes = merge((; color=color), attributes)
        band_color = get(band_attributes, :color, :black)
        band_attributes = merge((; color=band_color), band_attributes)
        lines_color = get(attributes, :color, :black)
        lines_attributes = merge((; color=lines_color), lines_attributes)
    end

    attributes = _extract_title!(ax, attributes)
    ax.xlabel = xlabel

    for (idx, (bootstrap, label)) in enumerate(zip(xs, labels))
        df = transform!(DataFrame(bootstrap.β), :coefname => ByRow(string) => :coefname)
        filter!(:coefname => in(_coefnames(bootstrap; show_intercept)), df)
        gdf = groupby(df, :coefname)
        dens = combine(gdf, :β => kde => :kde)

        if !ismissing(conf_level)
            coefplot!(ax, bootstrap;
                      conf_level,
                      vline_at_zero,
                      show_intercept,
                      show_legend=false,
                      color=Cycled(idx),
                      labels=[label],
                      scatter_attributes,
                      errorbars_attributes,
                      attributes...)
        end

        for (offset, row) in enumerate(reverse(eachrow(dens)))
            dd = 0.95 * row.kde.density ./ maximum(row.kde.density)
            lower = Point2f.(row.kde.x, offset)
            upper = Point2f.(row.kde.x, dd .+ offset)
            band!(ax, lower, upper;
                  color=Cycled(idx),
                  alpha=0.3,
                  attributes...,
                  band_attributes...,
                  label)
            lines!(ax, upper;
                   color=Cycled(idx),
                   attributes...,
                   lines_attributes...,
                   label)
        end
    end

    if ismissing(conf_level)
        nticks = _npreds(x; show_intercept)
        ax.yticks = (nticks:-1:1, cn)
        ylims!(ax, 0, nticks + 1)
        vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)
    end

    if show_legend
        legend_attributes = merge((; merge=true, unique=true), legend_attributes)
        axislegend(ax; legend_attributes...)
    end

    reset_limits!(ax)

    return ax
end
