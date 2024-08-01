"""
    ridgeplot(x::Union{MixedModel,MixedModelBootstrap}; kwargs...)::Figure
    ridgeplot!(fig::$(Indexable), x::Union{MixedModel,MixedModelBootstrap};
              kwargs...)
    ridgeplot!(ax::Axis, Union{MixedModel,MixedModelBootstrap};
              conf_level=0.95, vline_at_zero=true, show_intercept=true,
              scatter_attributes=(;),
              errorbars_attributes=(;),
              band_attributes=(;),
              lines_attributes=(;),
              attributes...)

Create a ridge plot for the bootstrap samples of the fixed effects.

Densities are normalized so that the maximum density is always 1.

The highest density interval corresponding to `conf_level` is marked with a bar at the bottom of each density.
Setting `conf_level=missing` removes the markings for the highest density interval.

`attributes` are passed onto [`coefplot`](@ref), `band!` and `lines!`.
`scatter_attributes` and `errorbars_attributes` are passed only onto [`coefplot`](@ref).
`band_attributes` and `lines_attributes` are passed only onto `band!` and
`linesbars!`, respectively.
(Starting with Makie 0.21, unsupported attributes for a
given plottype are no longer silently ignored, so it's necessary so separate out the
attributes that are only valid for a single plot type.)

The mutating methods return the original object.

!!! note
    Inestimable coefficients (coefficients removed by pivoting in the rank deficient case)
    are excluded.
"""
function ridgeplot(x::MixedModelBootstrap; show_intercept=true, kwargs...)
    # need to guarantee a min height of 200
    fig = Figure(; size=(640, max(200, 100 * _npreds(x; show_intercept))))
    return ridgeplot!(fig, x; show_intercept, kwargs...)
end

"""$(@doc ridgeplot)"""
function ridgeplot!(fig::Indexable, x::MixedModelBootstrap; kwargs...)
    ax = Axis(fig[1, 1])
    ridgeplot!(ax, x; kwargs...)
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
function ridgeplot!(ax::Axis, x::MixedModelBootstrap;
                    conf_level=0.95,
                    vline_at_zero=true,
                    show_intercept=true,
                    scatter_attributes=(;),
                    errorbars_attributes=(;),
                    band_attributes=(;),
                    lines_attributes=(;),
                    attributes...)
    xlabel = if !ismissing(conf_level)
        @sprintf "Normalized bootstrap density and %g%% confidence interval" (conf_level *
                                                                              100)
    else
        "Normalized bootstrap density"
    end

    if !ismissing(conf_level)
        coefplot!(ax, x; conf_level, vline_at_zero, show_intercept, color=:black,
                  scatter_attributes, errorbars_attributes, attributes...)
    end

    attributes = merge((; color=:black), attributes)
    band_attributes = merge((; color=(_color(attributes.color), 0.3)), band_attributes)

    ax.xlabel = xlabel

    df = transform!(DataFrame(x.β), :coefname => ByRow(string) => :coefname)
    filter!(:coefname => in(_coefnames(x; show_intercept)), df)
    group = :coefname
    densvar = :β
    gdf = groupby(df, group)

    dens = combine(gdf, densvar => kde => :kde)

    for (offset, row) in enumerate(reverse(eachrow(dens)))
        # multiply by 0.95 so that the ridges don't overlap
        dd = 0.95 * row.kde.density ./ maximum(row.kde.density)
        lower = Point2f.(row.kde.x, offset)
        upper = Point2f.(row.kde.x, dd .+ offset)
        band!(ax, lower, upper; attributes..., band_attributes...)
        lines!(ax, upper; attributes..., lines_attributes...)
    end

    # check conf_level so that we don't double print
    # if coefplot took care of it for us
    if ismissing(conf_level)
        cn = _coefnames(x; show_intercept)
        nticks = _npreds(x; show_intercept)
        ax.yticks = (nticks:-1:1, cn)
        ylims!(ax, 0, nticks + 1)
        vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)
    end

    reset_limits!(ax)

    return ax
end
