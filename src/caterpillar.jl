"""
    RanefInfo

Information on random effects conditional modes/means, variances, etc.

Used for creating caterpillar plots.

!!! note
    This functionality may be moved upstream into MixedModels.jl in the near future.
"""
struct RanefInfo{T<:AbstractFloat}
    cnames::Vector{String}
    levels::Vector
    ranef::Matrix{T}
    stddev::Matrix{T}
end

"""
    ranefinfo(m::LinearMixedModel)

Return a `NamedTuple{fnames(m), NTuple(k, RanefInfo)}` from model `m`
"""
function ranefinfo(m::LinearMixedModel{T}) where {T}
    fn = fnames(m)
    val = sizehint!(RanefInfo[], length(fn))
    for (re, eff, cv) in zip(m.reterms, ranef(m), condVar(m))
        push!(
            val,
            RanefInfo(
                re.cnames,
                re.levels,
                Matrix(adjoint(eff)),
                Matrix(adjoint(dropdims(sqrt.(mapslices(diag, cv, dims=1:2)), dims=2))),
            )
        )
    end
    NamedTuple{fn}((val...,))
end

"""
    caterpillar!(f::Figure, r::RanefInfo; orderby=1)

Add Axes of a caterpillar plot from `r` to `f`.

The order of the levels on the vertical axes is increasing `orderby` column
of `r.ranef`, usually the `(Intercept)` random effects.
Setting `orderby=nothing` will disable sorting, i.e. return the levels in the
order they are stored in.
"""
function caterpillar!(f::Figure, r::RanefInfo; orderby=1)
    rr = r.ranef
    if orderby === nothing
        ord = 1:size(rr, 1)
    else
        vv = view(rr, :, orderby)
        ord = sortperm(vv)
    end
    y = axes(rr, 1)
    cn = r.cnames
    axs = [Axis(f[1, j]) for j in axes(rr, 2)]
    linkyaxes!(axs...)
    for (j, ax) in enumerate(axs)
        xvals = view(rr, ord, j)
        scatter!(ax, xvals, y, color=(:red, 0.2))
        errorbars!(ax, xvals, y, view(r.stddev, ord, j), direction=:x)
        ax.xlabel = cn[j]
        ax.yticks = y
        j > 1 && hideydecorations!(ax, grid=false)
    end
    axs[1].yticks = (y, r.levels[ord])
    f
end

"""
    caterpillar(m::LinearMixedModel, gf::Symbol)

Returns a `Figure` of a "caterpillar plot" of the random-effects means and prediction intervals

A "caterpillar plot" is a horizontal error-bar plot of conditional means and standard deviations
of the random effects.
"""
function caterpillar(m::LinearMixedModel, gf::Symbol=first(fnames(m)))
    caterpillar!(Figure(resolution=(1000,800)), ranefinfo(m)[gf])
end
