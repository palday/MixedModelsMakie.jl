
@recipe(RidgePlot, x) do scene
    return Attributes(;
                      conf_level=0.95,
                      vline_at_zero=true,
                      show_intercept=true,
                      par=:β)
end

"""
    ridgeplot(x::MixedModelBootstrap;
              conf_level=0.95, vline_at_zero=true, show_intercept=true,
              attributes...)

Create a ridge plot for the bootstrap samples of the fixed effects.

Densities are normalized so that the maximum density is always 1.

The highest density interval correspoding to `conf_level` is marked with a bar at the bottom of each density.
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
    fig = Figure(; resolution=(640, max(200, 100 * _npreds(x; show_intercept))))
    ax = Axis(fig[1, 1])
    if !ismissing(conf_level)
        pl = coefplot!(ax, x; conf_level, vline_at_zero, show_intercept, attributes...)
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
    par = plot.par[]
    df = DataFrame(getproperty(boot, par))
    rename!(x -> replace(x, "column" => "coefname"), df)
    df = transform!(df, :coefname => ByRow(string) => :coefname)
    plot.show_intercept[] || filter!(:coefname => !=("(Intercept)"), df)
    group = :coefname
    # drop trailing s
    densvar = replace(string(par), "s" => "")
    # TODO: coefplot
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
