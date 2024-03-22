"""
    ridgeplot(x::Union{MixedModel,MixedModelBootstrap}; kwargs...)::Figure
    ridgeplot!(fig::$(Indexable), x::Union{MixedModel,MixedModelBootstrap}; 
              kwargs...)
    ridgeplot!(ax::Axis, Union{MixedModel,MixedModelBootstrap};
              conf_level=0.95, vline_at_zero=true, show_intercept=true, attributes...)

Create a ridge plot for the bootstrap samples of the fixed effects.

Densities are normalized so that the maximum density is always 1.

The highest density interval corresponding to `conf_level` is marked with a bar at the bottom of each density.
Setting `conf_level=missing` removes the markings for the highest density interval.

`attributes` are passed onto [`coefplot`](@ref), `band!` and `lines!`.

The mutating methods return the original object.
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
                    attributes...)
    xlabel = if !ismissing(conf_level)
        @sprintf "Normalized bootstrap density and %g%% confidence interval" (conf_level *
                                                                              100)
    else
        "Normalized bootstrap density"
    end

    if !ismissing(conf_level)
        coefplot!(ax, x; conf_level, vline_at_zero, show_intercept, color=:black,
                  attributes...)
    end

    attributes = merge((; xlabel, color=:black), attributes)
    band_attributes = merge(attributes, (; color=(_color(attributes.color), 0.3)))

    ax.xlabel = attributes.xlabel

    df = transform!(DataFrame(x.β), :coefname => ByRow(string) => :coefname)
    show_intercept || filter!(:coefname => !=("(Intercept)"), df)
    group = :coefname
    densvar = :β
    gdf = groupby(df, group)

    y = length(gdf):-1:1

    dens = combine(gdf, densvar => kde => :kde)

    for (offset, row) in enumerate(reverse(eachrow(dens)))
        # multiply by 0.95 so that the ridges don't overlap
        dd = 0.95 * row.kde.density ./ maximum(row.kde.density)
        lower = Point2f.(row.kde.x, offset)
        upper = Point2f.(row.kde.x, dd .+ offset)
        band!(ax, lower, upper; band_attributes...)
        lines!(ax, upper; attributes...)
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

# """
#     ridgeplot!(ax::Axis, df::AbstractDataFrame, densvar::Symbol, group::Symbol; normalize=false)
#     ridgeplot!(f::Union{Makie.FigureLike,Makie.GridLayout}, args...; pos=(1,1) kwargs...)
#     ridgeplot(args...; kwargs...)

# Create a "ridge plot".

# A ridge plot is stacked plot of densities for a given variable (`densvar`) grouped by a different variable (`group`). Because densities can very widely in scale, it is sometimes useful to `normalize` the densities so that each density has a maximum of 1.

# The non-mutating method creates a Figure before calling the method for Figure.
# The method for Figure places the ridge plot in the grid position specified by `pos`, default is (1,1).
# """
# function ridgeplot!(ax::Axis, df::AbstractDataFrame, densvar::Symbol, group::Symbol; sort_by_group=false, vline_at_zero=true, normalize=false)
#     # `normalize` makes it so that the max density is always 1
#     # `normalize` works on the density not the area/mass
#     gdf = groupby(df, group)
#     dens = combine(gdf, densvar => kde => :kde)
#     sort_by_group && sort!(dens, group)
#     spacing = normalize ? 1.0 :  0.9 * maximum(dens[!, :kde]) do val
#         return maximum(val.density)
#     end

#     nticks = length(gdf)

#     for (idx, row) in enumerate(eachrow(dens))
#         dd = normalize ? row.kde.density ./ maximum(row.kde.density) : row.kde.density

#         offset =  idx * spacing

#         lower = Node(Point2f.(row.kde.x, offset))
#         upper = Node(Point2f.(row.kde.x, dd .+ offset))
#         band!(ax, lower, upper; color=(:black, 0.3))
#         lines!(ax, upper; color=(:black, 1.0))
#     end

#     vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)

#     ax.yticks[] = (1:spacing:(nticks*spacing), string.(dens[!, group]))
#     ylims!(ax, 0, (nticks + 2) * spacing)
#     ax.xlabel[] = string(densvar)
#     ax.ylabel[] = string(group)

#     ax
# end

# function ridgeplot!(f::Union{Makie.FigureLike,Makie.GridLayout}, args...; pos=(1,1), kwargs...)
#     ridgeplot!(Axis(f[pos...]), args...; kwargs...)
#     return f
# end

# """
#     ridgeplot(args...; kwargs...)

# See [ridgeplot!](@ref).
# """
# ridgeplot(args...; kwargs...) = ridgeplot!(Figure(), args...; kwargs...)

# """
#     ridgeplot!(Union{Figure, Axis}, bstrp::MixedModelBootstrap, args...; show_intercept=true, kwargs...)
#     ridgeplot(bstrp::MixedModelBootstrap, args...; show_intercept=true, kwargs...)

# Convenience methods that call `DataFrame(bstrp.β)` and pass that onto the primary `ridgeplot[!]` methods.
# By default, the intercept is shown, but this can be disabled.
# """
# function ridgeplot!(axis::Axis, bstrp::MixedModelBootstrap, args...;
#                     show_intercept=true,  normalize=true, kwargs...)

#     df = DataFrame(bstrp.β)
#     show_intercept || filter!(:coefname => !=(Symbol("(Intercept)")), df)
#     ridgeplot!(axis, df, :β, :coefname, args...; sort_by_group=false, normalize, kwargs...)
# end
