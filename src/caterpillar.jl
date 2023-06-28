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
    ranefinfo(m::MixedModel)

Return a `NamedTuple{fnames(m), NTuple(k, RanefInfo)}` from model `m`
"""
function ranefinfo(m::MixedModel{T}) where {T}
    fn = fnames(m)
    val = sizehint!(RanefInfo[], length(fn))
    re = ranef(m)
    for grp in fn
        push!(val, ranefinfo(m, grp, re))
    end
    return NamedTuple{fn}((val...,))
end

"""
    ranefinfo(m::MixedModel, gf::Symbol)

Return a `RanefInfo` corresponding to the grouping variable `gf` model `m`.
"""
function ranefinfo(m::LinearMixedModel, gf::Symbol, re=ranef(m))
    idx = findfirst(==(gf), fnames(m))
    isnothing(idx) &&
        throw(ArgumentError("$gf is not the name of a grouping variable in the model"))

    # XXX replace ranef(m)[idx] with ranef(m, gf) when that becomes available upstream
    re, eff, cv = m.reterms[idx], re[idx], condVar(m, gf)
    return RanefInfo(re.cnames,
                     re.levels,
                     Matrix(adjoint(eff)),
                     Matrix(adjoint(dropdims(sqrt.(mapslices(diag, cv; dims=1:2)); dims=2))))
end

function ranefinfo(m::GeneralizedLinearMixedModel, args...; kwargs...)
    return ranefinfo(m.LMM, args...; kwargs...)
end

"""
    ranefinfotable(ri::RanefInfo)

Return the information in `ri` as a column table (`NamedTuple` of `Vector`s)

The columns are

- `name`: name of the random effect
- `level`: level of the grouping factor
- `cmode`: conditional mode of the random effect
- `cstddev`: conditional standard deviation of the random effect

"""
function ranefinfotable(ri::RanefInfo)
    cnames, levels = ri.cnames, ri.levels
    k = length(cnames)
    l = length(levels)
    return (;
            name=repeat(cnames; inner=l),
            level=repeat(levels; outer=k),
            cmode=vec(ri.ranef),
            cstddev=vec(ri.stddev))
end

"""
    caterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, r::RanefInfo;
                orderby=1, cols=nothing)
    caterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, m::LinearMixedModel,
                 gf::Symbol=first(fnames(m)); orderby=1, cols=nothing)

Add Axes of a caterpillar plot from `r` to `f`.

The order of the levels on the vertical axes is increasing `orderby` column
of `r.ranef`, usually the `(Intercept)` random effects.
Setting `orderby=nothing` will disable sorting, i.e. return the levels in the
order they are stored in.

The display can be restricted to a subset of random effects associated with a grouping variable by
specifying `cols`, either by indices or term names.

!!! note
    Even when not sorting the levels, they might have already been sorted during
    model matrix construction. If you want impose a particular ordering on the
    levels, then you must sort the relevant fields in the `RanefInfo` object before
    calling `caterpillar!`.

!!! note
    `orderby` is the ``n``th column of the columns specified by `cols`.
"""
function caterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, r::RanefInfo;
                      orderby=1, cols=nothing)
    cols = something(cols, axes(r.cnames, 1))
    cols = _cols_to_idx(r, cols)
    rr = view(r.ranef, :, cols)
    sd = view(r.stddev, :, cols)
    y = axes(rr, 1)
    ord = isnothing(orderby) ? y : sortperm(view(rr, :, orderby))
    cn = r.cnames
    axs = [Axis(f[1, j]) for j in axes(rr, 2)]
    linkyaxes!(axs...)
    for (j, ax) in enumerate(axs)
        xvals = view(rr, ord, j)
        scatter!(ax, xvals, y; color=(:red, 0.2))
        errorbars!(ax, xvals, y, 1.960 * view(sd, ord, j); direction=:x)
        ax.xlabel = cn[j]
        ax.yticks = y
        j > 1 && hideydecorations!(ax; grid=false)
    end
    axs[1].yticks = (y, string.(r.levels[ord]))
    return f
end

function caterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, m::MixedModel,
                      gf::Symbol=first(fnames(m)); kwargs...)
    return caterpillar!(f, ranefinfo(m, gf); kwargs...)
end

"""
    caterpillar(m::LinearMixedModel, gf::Symbol; orderby=1, cols=nothing)

Returns a `Figure` of a "caterpillar plot" of the random-effects means and prediction intervals

A "caterpillar plot" is a horizontal error-bar plot of conditional means and standard deviations
of the random effects.

The order of the levels on the vertical axes is increasing `orderby` column
of `r.ranef`, usually the `(Intercept)` random effects.
Setting `orderby=nothing` will disable sorting, i.e. return the levels in the
order they are stored in.

The display can be restricted to a subset of random effects associated with a grouping variable by
specifying `cols`, either by indices or term names.

!!! note
    `orderby` is the ``n``th column of the columns specified by `cols`.
"""
function caterpillar(m::MixedModel, gf::Symbol=first(fnames(m)); kwargs...)
    return caterpillar!(Figure(; resolution=(1000, 800)), m; kwargs...)
end

"""
    qqcaterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, r::RanefInfo; cols=nothing)
    qqcaterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, m::LinearMixedModel,
                   gf::Symbol=first(fnames(m)); cols=nothing)

Update the figure with a caterpillar plot with the vertical axis on the Normal() quantile scale.

The display can be restricted to a subset of random effects associated with a grouping variable by
specifying `cols`, either by indices or term names.

The order of the levels on the vertical axes is increasing `orderby` column
of `r.ranef`, usually the `(Intercept)` random effects.
Setting `orderby=nothing` will disable sorting, i.e. return the levels in the
order they are stored in.
"""
function qqcaterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, r::RanefInfo;
                        cols=nothing)
    cols = something(cols, axes(r.cnames, 1))
    cols = _cols_to_idx(r, cols)
    cn, rr = r.cnames, r.ranef
    y = zquantile.(ppoints(size(rr, 1)))
    axs = [Axis(f[1, j]) for j in axes(cols, 1)]
    linkyaxes!(axs...)
    for (j, k) in enumerate(cols)
        ax = axs[j]
        xvals = rr[:, k]
        ord = sortperm(xvals)
        xvals = xvals[ord]
        scatter!(ax, xvals, y; color=(:red, 0.2))
        errorbars!(ax, xvals, y, 1.960 * view(r.stddev, ord, k); direction=:x)
        ax.xlabel = string(cn[k])
        j > 1 && hideydecorations!(ax; grid=false)
    end
    return f
end

function qqcaterpillar!(f::Union{Makie.FigureLike,Makie.GridLayout}, m::MixedModel,
                        gf::Symbol=first(fnames(m)); kwargs...)
    return qqcaterpillar!(f, ranefinfo(m, gf); kwargs...)
end

"""
    qqcaterpillar(m::LinearMixedModel, gf::Symbol=first(fnames(m)); cols=nothing, orderby=1)

Returns a `Figure` of a "qq-caterpillar plot" of the random-effects means and prediction intervals.

The display can be restricted to a subset of random effects associated with a grouping variable by
specifying `cols`, either by indices or term names.
"""
function qqcaterpillar(m::MixedModel, gf::Symbol=first(fnames(m)); kwargs...)
    return qqcaterpillar!(Figure(; resolution=(1000, 800)), m, gf; kwargs...)
end

_cols_to_idx(r, cols) = cols
_cols_to_idx(r, cols::AbstractVector{<:Symbol}) = _cols_to_idx(r, string.(cols))
function _cols_to_idx(r, cols::Vector{<:AbstractString})
    return [i for (i, c) in enumerate(r.cnames) if c in cols]
end
