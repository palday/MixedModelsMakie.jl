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
            yrange=(minimum(yall), maximum(yall)))
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

"""
    facetregression(data, response, predictor, group; kwargs...)::Figure
    facetregression!(f::$(Indexable), data, response, predictor, group;
                     orderby::Symbol=:none,
                     rev::Bool=false,
                     layout::Union{Nothing,Tuple{Union{Nothing,Int},Union{Nothing,Int}}}=nothing,
                     bank45::Bool=true,
                     xlabel::AbstractString=string(predictor),
                     ylabel::AbstractString=string(response),
                     scattercolor=(:blue, 0.4),
                     linecolor=:red,
                     labelcolor=:black,
                     labelsize=12)

Create a Cleveland-trellis-style small-multiples display: one panel per level of
`group`, each showing a scatter of `response` against `predictor` plus that
group's own OLS regression line (via [`simplelinreg`](@ref)).

`response`, `predictor`, and `group` may be `Symbol` or `AbstractString` column
names into any Tables.jl-compatible `data`.

All panels share linked x/y axes (Cleveland-trellis convention), so slopes and
scatter remain visually comparable across panels.

# Keywords
- `orderby::Symbol=:none`: panel order — `:none` preserves the order groups are
  first encountered in `data` (not alphabetical), `:intercept`/`:slope` sort by
  the group's own OLS fit.
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
  per-group OLS slopes appears at 45° on screen ("banking", after Cleveland's
  trellis displays), easing visual comparison of slopes across panels. Full
  data is always shown — banking reshapes the panels, it never clips data.
  `bank45=false` leaves panels at their natural (roughly square) shape.
- `xlabel`/`ylabel`: shared axis titles spanning the whole grid (default to the
  `predictor`/`response` column names).
- `scattercolor`, `linecolor`: per-panel scatter and regression-line colors.
- `labelcolor`, `labelsize`: styling for the per-panel group-level corner label.

The mutating methods return the original object.
"""
function _facetregression_render!(f::Indexable, info::NamedTuple;
                                  orderby::Symbol=:none,
                                  rev::Bool=false,
                                  layout::Union{Nothing,Tuple{_MaybeInt,_MaybeInt}}=nothing,
                                  bank45::Bool=true,
                                  xlabel::AbstractString=string(info.predictor),
                                  ylabel::AbstractString=string(info.response),
                                  scattercolor=(:blue, 0.4),
                                  linecolor=:red,
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

function facetregression!(f::Indexable, data, response::Union{Symbol,AbstractString},
                          predictor::Union{Symbol,AbstractString},
                          group::Union{Symbol,AbstractString}; kwargs...)
    df = DataFrame(data)
    info = _facetregression_data(df, Symbol(response), Symbol(predictor), Symbol(group))
    return _facetregression_render!(f, info; kwargs...)
end

"""$(@doc facetregression!)"""
function facetregression(data, response, predictor, group; kwargs...)
    return facetregression!(Figure(; size=(1000, 800)), data, response, predictor, group;
                            kwargs...)
end

"""
    facetregressiontable(data, response, predictor, group; orderby::Symbol=:none,
                         rev::Bool=false)::DataFrame

Return the per-group OLS fits underlying [`facetregression!`](@ref) as a
`DataFrame`, one row per level of `group`, with columns:
- `group`: the group level
- `n`: number of observations in that group
- `intercept`, `slope`: the group's OLS fit of `response` on `predictor`
  (via [`simplelinreg`](@ref))

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
