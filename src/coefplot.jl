"""
    coefplot(x::MixedModel;
             conf_level=0.95, vline_at_zero=true, show_intercept=true, attributes...)
    coefplot(x::MixedModelBootstrap;
             conf_level=0.95, vline_at_zero=true, show_intercept=true, attributes...)

Create a coefficient plot of the fixed-effects and associated confidence intervals.

!!! note
    This functionality is implemented using [Makie recipes](https://makie.juliaplots.org/v0.15.0/recipes.html)
    and thus there are also additional auto-generated methods for `coefplot` and `coefplot!` that may be useful
    when constructing more complex figures.
"""
function coefplot(x::Union{MixedModel,MixedModelBootstrap};
                  conf_level=0.95,
                  vline_at_zero=true,
                  show_intercept=true,
                  attributes...)
    # need to guarantee a min height of 150
    fig = Figure(; resolution=(640, max(150, 75 * _npreds(x; show_intercept))))
    ax = Axis(fig[1, 1])
    pl = coefplot!(ax, x; conf_level, vline_at_zero, show_intercept, attributes...)
    return Makie.FigureAxisPlot(fig, ax, pl)
end

@recipe(CoefPlot, x) do scene
    return Attributes(; conf_level=0.95, vline_at_zero=true, show_intercept=true)
end

function Makie.plot!(ax::Axis, P::Type{<:CoefPlot}, allattrs::Makie.Attributes, x)
    plot = Makie.plot!(ax.scene, P, allattrs, x)

    if haskey(allattrs, :title)
        ax.title = allattrs.title[]
    end
    if haskey(allattrs, :xlabel)
        ax.xlabel = allattrs.xlabel[]
    else
        ax.xlabel = @sprintf "Estimate and %g%% confidence interval" (allattrs.conf_level[] *
                                                                      100)
    end
    if haskey(allattrs, :ylabel)
        ax.ylabel = allattrs.ylabel[]
    end
    reset_limits!(ax)
    show_intercept = allattrs.show_intercept[]
    cn = _coefnames(x; show_intercept)
    nticks = _npreds(x; show_intercept)
    ax.yticks = (nticks:-1:1, cn)
    ylims!(ax, 0, nticks + 1)
    allattrs.vline_at_zero[] && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)
    return plot
end

function Makie.plot!(plot::CoefPlot{<:Tuple{Union{MixedModel,MixedModelBootstrap}}})
    model_or_boot, conf_level = plot[1][], plot.conf_level[]
    ci = confint_table(model_or_boot, conf_level)
    plot.show_intercept[] || filter!(:coefname => !=("(Intercept)"), ci)
    y = nrow(ci):-1:1
    xvals = ci.estimate
    scatter!(plot, xvals, y)
    errorbars!(plot, xvals, y, xvals .- ci.lower, ci.upper .- xvals; direction=:x)

    return plot
end
