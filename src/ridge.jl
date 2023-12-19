
@recipe(RidgePlot, x) do scene
    return Attributes(;
                      conf_level=0.95,
                      vline_at_zero=true,
                      show_intercept=true)
end

"""
    ridgeplot(x::MixedModelBootstrap;
              conf_level=0.95, vline_at_zero=true, show_intercept=true,
              attributes...)

Create a ridge plot for the bootstrap samples of the fixed effects.

Densities are normalized so that the maximum density is always 1.

The highest density interval corresponding to `conf_level` is marked with a bar at the bottom of each density.
Setting `conf_level=missing` removes the markings for the highest density interval.

!!! note
    This functionality is implemented using [Makie recipes](https://makie.juliaplots.org/v0.15.0/recipes.html)
    and thus there are also additional auto-generated methods for `ridgeplot` and `ridgeplot!` that may be useful
    when constructing more complex figures.
"""
function ridgeplot(x::MixedModelBootstrap;
                   conf_level=0.95,
                   vline_at_zero=true,
                   show_intercept=true,
                   attributes...)
    # need to guarantee a min height of 200
    fig = Figure(; size=(640, max(200, 100 * _npreds(x; show_intercept))))
    ax = Axis(fig[1, 1])
    if !ismissing(conf_level)
        pl = coefplot!(ax, x; conf_level, vline_at_zero, show_intercept, color=:black, attributes...)
    end
    pl = ridgeplot!(ax, x; vline_at_zero, conf_level, show_intercept, attributes...)
    return Makie.FigureAxisPlot(fig, ax, pl)
end

function Makie.plot!(ax::Axis, P::Type{<:RidgePlot}, allattrs::Makie.Attributes, x)
    plot = Makie.plot!(ax.scene, P, allattrs, x)

    if haskey(allattrs, :title)
        ax.title = allattrs.title[]
    end
    if haskey(allattrs, :xlabel)
        ax.xlabel = allattrs.xlabel[]
    else
        lab = if !ismissing(allattrs.conf_level[])
            @sprintf "Normalized bootstrap density and %g%% confidence interval" (allattrs.conf_level[] *
                                                                                  100)
        else
            "Normalized bootstrap density"
        end
        ax.xlabel = lab
    end
    if haskey(allattrs, :ylabel)
        ax.ylabel = allattrs.ylabel[]
    end
    reset_limits!(ax)
    show_intercept = allattrs.show_intercept[]

    # check conf_level so that we don't double print
    # if coefplot took care of it for us
    if ismissing(allattrs.conf_level[])
        cn = _coefnames(x; show_intercept)
        nticks = _npreds(x; show_intercept)
        ax.yticks = (nticks:-1:1, cn)
        ylims!(ax, 0, nticks + 1)
        allattrs.vline_at_zero[] && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)
    end
    return plot
end

function Makie.plot!(plot::RidgePlot{<:Tuple{MixedModelBootstrap}})
    boot, conf_level = plot[1][], plot.conf_level[]
    df = transform!(DataFrame(boot.β), :coefname => ByRow(string) => :coefname)
    plot.show_intercept[] || filter!(:coefname => !=("(Intercept)"), df)
    group = :coefname
    densvar = :β
    gdf = groupby(df, group)

    y = length(gdf):-1:1

    dens = combine(gdf, densvar => kde => :kde)

    for (offset, row) in enumerate(reverse(eachrow(dens)))
        # multiply by 0.95 so that the ridges don't overlap
        dd = 0.95 * row.kde.density ./ maximum(row.kde.density)
        lower = Observable(Point2f.(row.kde.x, offset))
        upper = Observable(Point2f.(row.kde.x, dd .+ offset))
        band!(plot, lower, upper; color=(:black, 0.3))
        lines!(plot, upper; color=(:black, 1.0))
    end

    return plot
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
