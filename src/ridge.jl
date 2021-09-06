"""
    ridgeplot!(ax::Axis, df::AbstractDataFrame, densvar::Symbol, group::Symbol; normalize=false)
    ridgeplot!(f::Figure, args...; pos=(1,1) kwargs...)
    ridgeplot(args...; kwargs...)

Create a "ridge plot".

A ridge plot is stacked plot of densities for a given variable (`densvar`) grouped by a different variable (`group`). Because densities can very widely in scale, it is sometimes useful to `normalize` the densities so that each density has a maximum of 1.

The non-mutating method creates a Figure before calling the method for Figure.
The method for Figure places the ridge plot in the grid position specified by `pos`, default is (1,1).
"""
function ridgeplot!(ax::Axis, df::AbstractDataFrame, densvar::Symbol, group::Symbol; sort_by_group=false, vline_at_zero=true, normalize=false)
    # `normalize` makes it so that the max density is always 1
    # `normalize` works on the density not the area/mass
    gdf = groupby(df, group)
    dens = combine(gdf, densvar => kde => :kde)
    sort_by_group && sort!(dens, group)
    spacing = normalize ? 1.0 :  0.9 * maximum(dens[!, :kde]) do val
        return maximum(val.density)
    end

    nticks = length(gdf)

    for (idx, row) in enumerate(eachrow(dens))
        dd = normalize ? row.kde.density ./ maximum(row.kde.density) : row.kde.density

        offset =  idx * spacing

        lower = Node(Point2f.(row.kde.x, offset))
        upper = Node(Point2f.(row.kde.x, dd .+ offset))
        band!(ax, lower, upper; color=(:black, 0.3))
        lines!(ax, upper; color=(:black, 1.0))
    end

    vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)

    ax.yticks[] = (1:spacing:(nticks*spacing), string.(dens[!, group]))
    ylims!(ax, 0, (nticks + 2) * spacing)
    ax.xlabel[] = string(densvar)
    ax.ylabel[] = string(group)

    ax
end

function ridgeplot!(f::Figure, args...; pos=(1,1), kwargs...)
    ridgeplot!(Axis(f[pos...]), args...; kwargs...)
    return f
end

"""
    ridgeplot(args...; kwargs...)

See [ridgeplot!](@ref).
"""
ridgeplot(args...; kwargs...) = ridgeplot!(Figure(), args...; kwargs...)

"""
    ridgeplot!(Union{Figure, Axis}, bstrp::MixedModelBootstrap, args...; show_intercept=true, kwargs...)
    ridgeplot(bstrp::MixedModelBootstrap, args...; show_intercept=true, kwargs...)

Convenience methods that call `DataFrame(bstrp.β)` and pass that onto the primary `ridgeplot[!]` methods.
By default, the intercept is shown, but this can be disabled.
"""
function ridgeplot!(axis::Axis, bstrp::MixedModelBootstrap, args...;
                    show_intercept=true, kwargs...)

    df = DataFrame(bstrp.β)
    hide_intercept && filter!(:coefname => !=(Symbol("(Intercept)")), df)
    ridgeplot!(axis, :β, :coefname; sort_by_group=false)
end
