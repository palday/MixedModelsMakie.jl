function _facetregression_data(df::DataFrame, response::Symbol, predictor::Symbol,
                               group::Symbol)
    xall = Float64.(df[!, predictor])
    yall = Float64.(df[!, response])
    gcol = df[!, group]
    labels = unique(gcol)
    x = Vector{Float64}[]
    y = Vector{Float64}[]
    n = Int[]
    intercept = Float64[]
    slope = Float64[]
    for lab in labels
        idx = findall(==(lab), gcol)
        xg, yg = xall[idx], yall[idx]
        a, b = simplelinreg(xg, yg)
        push!(x, xg)
        push!(y, yg)
        push!(n, length(idx))
        push!(intercept, a)
        push!(slope, b)
    end
    return (; group, predictor, response, labels, x, y, n, intercept, slope,
            xrange=(minimum(xall), maximum(xall)),
            yrange=(minimum(yall), maximum(yall)),
            fixef=nothing, shrunken_intercept=nothing, shrunken_slope=nothing)
end

"""
    _facetregression_data(m::LinearMixedModel, predictor, group)

Extract the same per-group data as the table-based method, sourced entirely
from the model's own matrices (`m.X`, `m.y`, `m.reterms`) rather than a raw
data table -- see `ranefinfo` (caterpillar.jl) and the co-occurrence
extraction in nesting.jl for the precedent of using only `ReMat.refs`/`.levels`
this way. Also populates `fixef` (the population intercept/slope) and, per
group, `shrunken_intercept`/`shrunken_slope` (fixed effects + that group's
conditional modes from `ranef(m)`).
"""
function _facetregression_data(m::LinearMixedModel, predictor::Union{Symbol,AbstractString},
                               group::Union{Symbol,AbstractString})
    group = Symbol(group)
    gidx = findfirst(==(group), fnames(m))
    gidx === nothing &&
        throw(ArgumentError("$group is not the name of a grouping variable in the model; available: $(fnames(m))"))
    re = m.reterms[gidx]
    labels = collect(re.levels)
    refs = re.refs
    nlevels = length(labels)

    cn = coefnames(m)
    predictor = string(predictor)
    pidx = findfirst(==(predictor), cn)
    pidx === nothing &&
        throw(ArgumentError("$predictor is not a fixed-effect term in the model; available: $cn"))
    iidx = findfirst(==("(Intercept)"), cn)
    iidx === nothing &&
        throw(ArgumentError("model has no fixed intercept; facetregression requires one"))

    xall = view(m.X, :, pidx)
    yall = m.y

    x = Vector{Float64}[]
    y = Vector{Float64}[]
    n = Int[]
    intercept = Float64[]
    slope = Float64[]
    for ℓ in 1:nlevels
        obs_idx = findall(==(ℓ), refs)
        xg, yg = xall[obs_idx], yall[obs_idx]
        a, b = simplelinreg(xg, yg)
        push!(x, xg)
        push!(y, yg)
        push!(n, length(obs_idx))
        push!(intercept, a)
        push!(slope, b)
    end

    fe = fixef(m)
    a_pop, b_pop = fe[iidx], fe[pidx]

    remat = ranef(m)[gidx]                                # (ncoef, nlevels)
    rint = findfirst(==("(Intercept)"), re.cnames)
    rslope = findfirst(==(predictor), re.cnames)
    shrunken_intercept = a_pop .+ (rint === nothing ? zeros(nlevels) : remat[rint, :])
    shrunken_slope = b_pop .+ (rslope === nothing ? zeros(nlevels) : remat[rslope, :])

    response = try
        Symbol(m.formula.lhs.sym)
    catch
        :response
    end

    return (; group, predictor=Symbol(predictor), response, labels, x, y, n, intercept,
            slope, xrange=(minimum(xall), maximum(xall)),
            yrange=(minimum(yall), maximum(yall)), fixef=(a_pop, b_pop),
            shrunken_intercept, shrunken_slope)
end

"""
    _facetregression_default_predictor(m::LinearMixedModel)

Return the model's only non-intercept fixed-effect term name, for the
`facetregression` methods that let `predictor` be omitted. Throws an
`ArgumentError` if the model doesn't have exactly one such term.
"""
function _facetregression_default_predictor(m::LinearMixedModel)
    candidates = filter(!=("(Intercept)"), coefnames(m))
    length(candidates) == 1 ||
        throw(ArgumentError("facetregression could not pick a default predictor: model has $(length(candidates)) non-intercept fixed effects ($candidates); pass `predictor` explicitly"))
    return only(candidates)
end

function _facetregression_order(info::NamedTuple, orderby::Symbol; rev::Bool=false)
    perm = if orderby === :none
        collect(eachindex(info.labels))
    elseif orderby === :intercept
        sortperm(info.intercept)
    elseif orderby === :slope
        sortperm(info.slope)
    else
        throw(ArgumentError("orderby must be :none, :intercept, or :slope, got :$orderby"))
    end
    return rev ? reverse(perm) : perm
end

const _MaybeInt = Union{Nothing,Int}

function _facetregression_layout(ngroups::Integer,
                                 layout::Union{Nothing,Tuple{_MaybeInt,_MaybeInt}})
    if layout === nothing
        ncol = ceil(Int, sqrt(ngroups))
        nrow = ceil(Int, ngroups / ncol)
        return (nrow, ncol)
    end
    nrow, ncol = layout
    if nrow === nothing && ncol === nothing
        ncol = ceil(Int, sqrt(ngroups))
        nrow = ceil(Int, ngroups / ncol)
    elseif nrow === nothing
        nrow = ceil(Int, ngroups / ncol)
    elseif ncol === nothing
        ncol = ceil(Int, ngroups / nrow)
    end
    return (nrow, ncol)
end

# `clevelandaxes!` only returns `f`, not the `GridLayout` it populated, and
# `f.layout` (when it exists) is the *parent* grid a GridPosition/
# GridSubposition lives in, not whatever got nested inside that cell -- so the
# actual layout holding the new axes must be found by walking the content tree.
_facetregression_gridlayout(f::Figure) = f.layout
_facetregression_gridlayout(f::GridLayout) = f
function _facetregression_gridlayout(f::Union{GridPosition,GridSubposition})
    idx = findfirst(c -> c isa GridLayout, Makie.contents(f))
    idx === nothing &&
        error("facetregression: expected clevelandaxes! to have created a nested GridLayout")
    return Makie.contents(f)[idx]
end

function _facetregression_render!(f::Indexable, info::NamedTuple;
                                  orderby::Symbol=:none,
                                  rev::Bool=false,
                                  layout::Union{Nothing,Tuple{_MaybeInt,_MaybeInt}}=nothing,
                                  bank45::Bool=true,
                                  xlabel::AbstractString=string(info.predictor),
                                  ylabel::AbstractString=string(info.response),
                                  scattercolor=(:blue, 0.4),
                                  linecolor=:red,
                                  show_fixef::Bool=true,
                                  fixefcolor=:black,
                                  fixeflinestyle=:dash,
                                  show_shrunken::Bool=true,
                                  shrunkencolor=:green,
                                  labelcolor=:black,
                                  labelsize=12)
    perm = _facetregression_order(info, orderby; rev)
    labels = info.labels[perm]
    lo = _facetregression_layout(length(labels), layout)
    nrow, ncol = lo

    npanel = nrow * ncol
    if npanel < length(labels)
        ndropped = length(labels) - npanel
        @warn "facetregression: layout=($nrow, $ncol) only fits $npanel of $(length(labels)) groups; $ndropped trailing group(s) will not be displayed"
        perm = perm[1:npanel]
        labels = labels[1:npanel]
    end

    clevelandaxes!(f, string.(labels), lo)
    gl = _facetregression_gridlayout(f)
    axs = filter(x -> x isa Axis, Makie.contents(gl))

    xmin, xmax = info.xrange
    dx_range = xmax - xmin
    xpad = 0.05 * dx_range
    xlo, xhi = xmin - xpad, xmax + xpad

    ymin, ymax = info.yrange
    dy_range = ymax - ymin
    ypad = 0.05 * dy_range
    ylo, yhi = ymin - ypad, ymax + ypad

    mean_slope = sum(info.slope) / length(info.slope)
    target_aspect = if bank45 && !iszero(mean_slope) && !iszero(dx_range) &&
                       !iszero(dy_range)
        abs(mean_slope) * dx_range / dy_range
    else
        bank45 &&
            @warn "facetregression: cannot bank to 45° (mean slope or data range is zero); using natural aspect"
        nothing
    end

    for (k, idx) in enumerate(perm)
        ax = axs[k]
        scatter!(ax, info.x[idx], info.y[idx]; color=scattercolor)
        xlims!(ax, xlo, xhi)
        ylims!(ax, ylo, yhi)
        ablines!(ax, info.intercept[idx], info.slope[idx]; color=linecolor)
        if show_fixef && info.fixef !== nothing
            a_pop, b_pop = info.fixef
            ablines!(ax, a_pop, b_pop; color=fixefcolor, linestyle=fixeflinestyle)
        end
        if show_shrunken && info.shrunken_intercept !== nothing
            ablines!(ax, info.shrunken_intercept[idx], info.shrunken_slope[idx];
                    color=shrunkencolor)
        end
        text!(ax, xlo + 0.05 * (xhi - xlo), yhi - 0.05 * (yhi - ylo);
             text=string(labels[k]), color=labelcolor,
             fontsize=labelsize, align=(:left, :top))
    end
    linkaxes!(axs...)

    # Shape the grid cells themselves to the banked aspect (rather than
    # constraining each Axis's own `aspect`, which only pads whitespace
    # inside an already-fixed-size cell and reopens the gaps `clevelandaxes!`
    # zeroed out).
    if target_aspect !== nothing
        for j in 1:ncol
            colsize!(gl, j, Aspect(1, target_aspect))
        end
    end

    # Explicit ranges (rather than `end`/`:`) so this works whether `f` is a
    # fresh Figure/GridLayout or a single fixed GridPosition/GridSubposition
    # cell, which don't support "current extent" indexing.
    Label(f[nrow + 1, 1:ncol], xlabel; tellwidth=false, tellheight=true)
    Label(f[1:nrow, 0], ylabel; tellwidth=true, tellheight=false, rotation=pi / 2)
    return f
end

"""
    facetregression(data, response, predictor, group; kwargs...)::Figure
    facetregression!(f::$(Indexable), data, response, predictor, group; kwargs...)
    facetregression(m::LinearMixedModel, predictor, group=first(fnames(m)); kwargs...)::Figure
    facetregression!(f::$(Indexable), m::LinearMixedModel, predictor,
                     group=first(fnames(m)); kwargs...)
    facetregression(m::LinearMixedModel; group=first(fnames(m)), kwargs...)::Figure
    facetregression!(f::$(Indexable), m::LinearMixedModel;
                     group=first(fnames(m)), kwargs...)

Create a Cleveland-trellis-style small-multiples display: one panel per level of
`group`, each showing a scatter of `response` against `predictor` plus that
group's own OLS regression line (via [`simplelinreg`](@ref)).

The table-based methods take `response`, `predictor`, and `group` as `Symbol` or
`AbstractString` column names into any Tables.jl-compatible `data`.

The `LinearMixedModel`-based methods instead extract everything from the model's
own matrices — no data table is needed. `predictor` must be the name of a
continuous fixed-effect term (e.g. `:days`); `group` names a grouping factor and
defaults to the first one (`first(fnames(m))`). Two extra reference lines are
drawn in every panel: the population-level fixed-effects fit (`show_fixef`,
the same line in every panel) and that group's own shrunken/BLUP fit
(`show_shrunken`, fixed effects plus that group's conditional modes from
`ranef(m)` — equal to the population fit if the model has no random slope on
`predictor` for `group`).

`predictor` may be omitted entirely (note `group` then becomes keyword-only, to
avoid ambiguity with the `predictor` positional argument), in which case it
defaults to the model's only non-intercept fixed-effect term — this throws an
`ArgumentError` if the model has more than one such term.

All panels share linked x/y axes (Cleveland-trellis convention), so slopes and
scatter remain visually comparable across panels.

# Keywords
- `orderby::Symbol=:none`: panel order — `:none` preserves the order groups are
  first encountered (not alphabetical), `:intercept`/`:slope` sort by the
  group's own within-unit OLS fit.
- `rev::Bool=false`: reverse the panel order produced by `orderby` (applied
  after sorting, so it also reverses `:none`'s first-encountered order).
- `layout::Union{Nothing,Tuple{Union{Nothing,Int},Union{Nothing,Int}}}=nothing`:
  `(nrow, ncol)` grid shape. `nothing` auto-computes a roughly square grid.
  Either element may be `nothing` to auto-compute just that one from the
  other (e.g. `layout=(2, nothing)` fixes 2 rows and picks enough columns to
  fit every group). If both are given and their product is smaller than the
  number of groups, the grid is used as given (never enlarged) and a warning
  reports how many trailing groups are left undisplayed.
- `bank45::Bool=true`: shape the grid columns (via `colsize!` with a
  `GridLayoutBase.Aspect` size) so that a line with the *mean* of the
  per-group within-unit OLS slopes appears at 45° on screen ("banking", after
  Cleveland's trellis displays), easing visual comparison of slopes across
  panels. Full data is always shown — banking reshapes the panels, it never
  clips data. `bank45=false` leaves panels at their natural (roughly square)
  shape.
- `xlabel`/`ylabel`: shared axis titles spanning the whole grid (default to the
  `predictor`/`response` names).
- `scattercolor`, `linecolor`: per-panel scatter and within-unit regression-line
  colors.
- `show_fixef::Bool=true`, `fixefcolor=:black`, `fixeflinestyle=:dash`: whether
  to draw the population fixed-effects line (`LinearMixedModel` source only).
- `show_shrunken::Bool=true`, `shrunkencolor=:green`: whether to draw the
  per-group shrunken/BLUP line (`LinearMixedModel` source only).
- `labelcolor`, `labelsize`: styling for the per-panel group-level corner label.

The mutating methods return the original object.
"""
function facetregression!(f::Indexable, data, response::Union{Symbol,AbstractString},
                          predictor::Union{Symbol,AbstractString},
                          group::Union{Symbol,AbstractString}; kwargs...)
    df = DataFrame(data)
    info = _facetregression_data(df, Symbol(response), Symbol(predictor), Symbol(group))
    return _facetregression_render!(f, info; kwargs...)
end

function facetregression!(f::Indexable, m::LinearMixedModel,
                          predictor::Union{Symbol,AbstractString},
                          group::Union{Symbol,AbstractString}=first(fnames(m)); kwargs...)
    info = _facetregression_data(m, predictor, group)
    return _facetregression_render!(f, info; kwargs...)
end

function facetregression!(f::Indexable, m::LinearMixedModel;
                          group::Union{Symbol,AbstractString}=first(fnames(m)), kwargs...)
    predictor = _facetregression_default_predictor(m)
    return facetregression!(f, m, predictor, group; kwargs...)
end

const _facetregression_doc = "$(@doc facetregression!(::Indexable, ::Any, ::Union{Symbol,AbstractString}, ::Union{Symbol,AbstractString}, ::Union{Symbol,AbstractString}))"

"""$(_facetregression_doc)"""
function facetregression(data, response, predictor, group; kwargs...)
    return facetregression!(Figure(; size=(1000, 800)), data, response, predictor, group;
                            kwargs...)
end

function facetregression(m::LinearMixedModel, predictor::Union{Symbol,AbstractString},
                         group::Union{Symbol,AbstractString}=first(fnames(m)); kwargs...)
    return facetregression!(Figure(; size=(1000, 800)), m, predictor, group; kwargs...)
end

function facetregression(m::LinearMixedModel;
                         group::Union{Symbol,AbstractString}=first(fnames(m)), kwargs...)
    return facetregression!(Figure(; size=(1000, 800)), m; group, kwargs...)
end

"""
    facetregressiontable(data, response, predictor, group; orderby::Symbol=:none)::DataFrame
    facetregressiontable(m::LinearMixedModel, predictor, group=first(fnames(m));
                         orderby::Symbol=:none)::DataFrame
    facetregressiontable(m::LinearMixedModel; group=first(fnames(m)),
                         orderby::Symbol=:none)::DataFrame

Return the per-group OLS fits underlying [`facetregression!`](@ref) as a
`DataFrame`, one row per level of `group`, with columns:
- `group`: the group level
- `n`: number of observations in that group
- `intercept`, `slope`: the group's within-unit OLS fit of `response` on `predictor`
  (via [`simplelinreg`](@ref))

As with [`facetregression!`](@ref), `predictor` may be omitted for the
`LinearMixedModel` methods (`group` then becomes keyword-only), defaulting to
the model's only non-intercept fixed-effect term.

The `LinearMixedModel` method additionally has:
- `fixef_intercept`, `fixef_slope`: the population fixed-effects fit (same for
  every row)
- `shrunken_intercept`, `shrunken_slope`: that group's shrunken/BLUP fit (fixed
  effects plus that group's conditional modes from `ranef(m)`)

`orderby` and `rev` have the same meaning as in [`facetregression!`](@ref).
"""
function facetregressiontable(data, response::Union{Symbol,AbstractString},
                              predictor::Union{Symbol,AbstractString},
                              group::Union{Symbol,AbstractString};
                              orderby::Symbol=:none, rev::Bool=false)
    df = DataFrame(data)
    info = _facetregression_data(df, Symbol(response), Symbol(predictor), Symbol(group))
    perm = _facetregression_order(info, orderby; rev)
    return DataFrame(; group=info.labels[perm], n=info.n[perm],
                     intercept=info.intercept[perm], slope=info.slope[perm])
end

function facetregressiontable(m::LinearMixedModel, predictor::Union{Symbol,AbstractString},
                              group::Union{Symbol,AbstractString}=first(fnames(m));
                              orderby::Symbol=:none, rev::Bool=false)
    info = _facetregression_data(m, predictor, group)
    perm = _facetregression_order(info, orderby; rev)
    a_pop, b_pop = info.fixef
    return DataFrame(; group=info.labels[perm], n=info.n[perm],
                     intercept=info.intercept[perm], slope=info.slope[perm],
                     fixef_intercept=fill(a_pop, length(perm)),
                     fixef_slope=fill(b_pop, length(perm)),
                     shrunken_intercept=info.shrunken_intercept[perm],
                     shrunken_slope=info.shrunken_slope[perm])
end

function facetregressiontable(m::LinearMixedModel;
                              group::Union{Symbol,AbstractString}=first(fnames(m)),
                              orderby::Symbol=:none, rev::Bool=false)
    predictor = _facetregression_default_predictor(m)
    return facetregressiontable(m, predictor, group; orderby, rev)
end
