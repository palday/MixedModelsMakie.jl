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
    fig = Figure(; size=(640, max(150, 75 * _npreds(x; show_intercept))))
    ax = Axis(fig[1, 1])
    coefplot!(ax, x; conf_level, vline_at_zero, show_intercept, attributes...)
    return Makie.FigureAxis(fig, ax)
end


# Indexable
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

    attributes = merge((;xlabel), attributes)

    scatter!(ax, xvals, y; attributes...)
    errorbars!(ax, xvals, y, xvals .- ci.lower, ci.upper .- xvals; direction=:x)
    vline_at_zero && vlines!(ax, 0; color=(:black, 0.75), linestyle=:dash)

    reset_limits!(ax)
    cn = _coefnames(x; show_intercept)
    nticks = _npreds(x; show_intercept)
    ax.yticks = (nticks:-1:1, cn)
    ylims!(ax, 0, nticks + 1)
    return ax
end
