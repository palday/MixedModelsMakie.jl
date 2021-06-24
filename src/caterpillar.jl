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

!!! note
    Even when not sorting the levels, they might have already been sorted during
    model matrix construction. If you want impose a particular ordering on the
    levels, then you must sort the relevant fields in the `RanefInfo` object before
    calling `caterpillar!`.
"""
function caterpillar!(f::Figure, r::RanefInfo; orderby=1)
    rr = r.ranef
    y = axes(rr, 1)
    ord = isnothing(orderby) ? y : sortperm(view(rr, :, orderby))
    cn = r.cnames
    axs = [Axis(f[1, j]) for j in axes(rr, 2)]
    linkyaxes!(axs...)
    for (j, ax) in enumerate(axs)
        xvals = view(rr, ord, j)
        scatter!(ax, xvals, y, color=(:red, 0.2))
        errorbars!(ax, xvals, y, 1.960 * view(r.stddev, ord, j), direction=:x)
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

"""
    qqcaterpillar!(f::Figure, r::RanefInfo; cols=axes(r.cnames, 1))

Update the figure with a caterpillar plot with the vertical axis on the Normal() quantile scale.

The display can be restricted to a subset of random effects associated with a grouping variable by specifying `cols`, either by indices or term names.
"""
function qqcaterpillar!(f::Figure, r::RanefInfo; cols=axes(r.cnames, 1))
    cols = _cols_to_idx(r, cols)
    cn, rr = r.cnames, r.ranef
    y = zquantile.(ppoints(size(rr, 1)))
    axs = [Axis(f[1,j]) for j in axes(cols, 1)]
    linkyaxes!(axs...)
    for (j, k) in enumerate(cols)
        ax = axs[j]
        xvals = rr[:, k]
        ord = sortperm(xvals)
        xvals = xvals[ord]
        scatter!(ax, xvals, y, color=(:red, 0.2))
        errorbars!(ax, xvals, y, 1.960 * view(r.stddev, ord, k); direction=:x)
        ax.xlabel = string(cn[k])
        j > 1 && hideydecorations!(ax, grid=false)
    end
    f
end

_cols_to_idx(r, cols) = cols
_cols_to_idx(r, cols::AbstractVector{<:Symbol}) = _cols_to_idx(r, string.(cols))
_cols_to_idx(r, cols::Vector{<:AbstractString}) = [i for (i,c) in enumerate(r.cnames) if c in cols]

"""
    qqcaterpillar(m::LinearMixedModel, gf::Symbol=first(fnames(m)); cols=nothing)

Returns a `Figure` of a "qq-caterpillar plot" of the random-effects means and prediction intervals.

The display can be restricted to a subset of random effects associated with a grouping variable by specifying `cols`, either by indices or term names.
"""
function qqcaterpillar(m::LinearMixedModel, gf::Symbol=first(fnames(m)); cols=nothing)
    reinfo = ranefinfo(m)[gf]
    cols = something(cols, axes(reinfo.cnames, 1))
    qqcaterpillar!(Figure(resolution=(1000, 800)), reinfo; cols=cols)
end
