"""
    _distinct_reterms(m::MixedModel) → Vector{Int}

Indices into `m.reterms` (equivalently `fnames(m)`) of the first occurrence of
each distinct grouping factor name, in order of first appearance.

A single grouping factor may back multiple random-effects terms, e.g.
`(1|subj) + (0+cond|subj)`. The nesting/crossing structure depends only on the
factor's level assignment, so duplicates are collapsed to one representative
term.
"""
function _distinct_reterms(m::MixedModel)
    fn = fnames(m)
    idxs = Int[]
    seen = Set{Symbol}()
    for (i, f) in enumerate(fn)
        f in seen && continue
        push!(seen, f)
        push!(idxs, i)
    end
    return idxs
end

"""
    _cooccurrence(refs_a, n_a, refs_b, n_b) → Matrix{Int}

Contingency table of co-occurrence counts between two grouping factors' level
assignments (one row per observation). `refs_a`/`refs_b` are integer group
indices, as stored in `ReMat.refs`; `n_a`/`n_b` are the number of levels of
each factor. The returned matrix has size `(n_a, n_b)`.
"""
function _cooccurrence(refs_a::AbstractVector{<:Integer}, n_a::Integer,
                       refs_b::AbstractVector{<:Integer}, n_b::Integer)
    tab = zeros(Int, n_a, n_b)
    for (a, b) in zip(refs_a, refs_b)
        tab[a, b] += 1
    end
    return tab
end

"""
    _nesting_relationship(tab::AbstractMatrix{<:Integer}) → NamedTuple

Classify the pairwise relationship between two grouping factors A (rows) and
B (columns) from their co-occurrence contingency table `tab`.

Returns `(; rel, density)` where `rel` is one of:
- `:identical`         — A and B partition observations identically (e.g. two
                          names for the same grouping factor)
- `:b_nested_in_a`     — every level of B occurs with exactly one level of A
- `:a_nested_in_b`     — every level of A occurs with exactly one level of B
- `:complete_crossing` — every combination of an A-level and a B-level is observed
- `:partial_crossing`  — some, but not all, combinations of A- and B-levels co-occur

`density` is the fraction of the `n_a × n_b` grid of combinations that is
actually observed (always `1.0` for `:complete_crossing`).

A dense, fully-crossed table can never also satisfy either nesting condition
once either factor has more than one level: nesting requires every row (or
column) to have a *single* nonzero entry, while complete crossing requires
*every* entry to be nonzero. The two conditions coincide only in the trivial
case of a permutation table, which is exactly the `:identical` case.
"""
function _nesting_relationship(tab::AbstractMatrix{<:Integer})
    n_a, n_b = size(tab)
    b_in_a = all(count(!iszero, view(tab, :, b)) == 1 for b in axes(tab, 2))
    a_in_b = all(count(!iszero, view(tab, a, :)) == 1 for a in axes(tab, 1))
    n_pairs = count(!iszero, tab)
    density = n_pairs / (n_a * n_b)
    rel = if a_in_b && b_in_a
        :identical
    elseif b_in_a
        :b_nested_in_a
    elseif a_in_b
        :a_nested_in_b
    elseif n_pairs == n_a * n_b
        :complete_crossing
    else
        :partial_crossing
    end
    return (; rel, density)
end

"""
    _block_order(tab::AbstractMatrix{<:Integer}) → (row_order, col_order)

Permutations of the rows and columns of `tab` that group each level with its
most common partner level, so that (near-)nested structure appears as a
block-diagonal pattern when `tab` is drawn as a heatmap. Purely a display aid;
it does not affect `_nesting_relationship`'s classification.
"""
function _block_order(tab::AbstractMatrix{<:Integer})
    n_a, n_b = size(tab)
    row_key = [argmax(view(tab, a, :)) for a in 1:n_a]
    col_key = [argmax(view(tab, :, b)) for b in 1:n_b]
    return sortperm(row_key), sortperm(col_key)
end

"""
    _nesting_data(m::MixedModel; gfs=nothing)

Build pairwise co-occurrence tables and relationship classifications for the
distinct grouping factors in `m`.

`gfs`, if given, restricts and orders the grouping factors to show (by name or
by index into `fnames(m)`); otherwise every distinct grouping factor is used,
in the order it first appears.

Returns a `NamedTuple` with fields `names` (`Vector{String}`), `nlevels`
(levels per factor), `levels` (level labels per factor, as `Vector{String}`),
`tables` (`Dict{Tuple{Int,Int},Matrix{Int}}`, keyed `(i, j)` with `i > j`,
rows indexed by factor `i` and columns by factor `j`), and `relationships`
(`Dict{Tuple{Int,Int},NamedTuple}`, same keys, as returned by
[`_nesting_relationship`](@ref)).
"""
function _nesting_data(m::MixedModel;
                       gfs::Union{Nothing,AbstractVector}=nothing)
    idxs = _distinct_reterms(m)
    if gfs !== nothing
        fn = fnames(m)
        wanted = [g isa Integer ? idxs[g] : findfirst(==(Symbol(g)), fn) for g in gfs]
        any(isnothing, wanted) &&
            throw(ArgumentError("gfs=$gfs contains a name that is not a grouping " *
                                "factor in the model; available names are $(unique(fn))"))
        idxs = wanted
    end
    length(idxs) < 2 &&
        throw(ArgumentError("Model has fewer than two distinct grouping factors " *
                            "to compare; nesting/crossing structure requires at least two."))

    names = string.(fnames(m)[idxs])
    reterms = m.reterms[idxs]
    levels = [string.(r.levels) for r in reterms]
    nlevels = length.(levels)

    tables = Dict{Tuple{Int,Int},Matrix{Int}}()
    relationships = Dict{Tuple{Int,Int},NamedTuple}()
    for i in 2:length(idxs), j in 1:(i - 1)
        tab = _cooccurrence(reterms[i].refs, nlevels[i], reterms[j].refs, nlevels[j])
        tables[(i, j)] = tab
        relationships[(i, j)] = _nesting_relationship(tab)
    end

    return (; names, nlevels, levels, tables, relationships)
end

"""
    _nesting_incidence_table(info::NamedTuple)

Convert the `info` NamedTuple (as returned by [`_nesting_data`](@ref)) into a
`DataFrame` in long format: one row per combination of a level of factor `A`
and a level of factor `B`, for every pair of grouping factors. Includes rows
with `count == 0` for combinations that never co-occur, since those absences
are exactly what reveal nesting or partial crossing.
"""
function _nesting_incidence_table(info::NamedTuple)
    factor_a = String[]
    level_a = String[]
    factor_b = String[]
    level_b = String[]
    count = Int[]

    for (i, j) in sort!(collect(keys(info.tables)))
        tab = info.tables[(i, j)]
        la, lb = info.levels[i], info.levels[j]
        for a in axes(tab, 1), b in axes(tab, 2)
            push!(factor_a, info.names[i])
            push!(level_a, la[a])
            push!(factor_b, info.names[j])
            push!(level_b, lb[b])
            push!(count, tab[a, b])
        end
    end

    return DataFrame(; factor_a, level_a, factor_b, level_b, count)
end

"""
    nestingtable(m::MixedModel; gfs::Union{Nothing,AbstractVector}=nothing)

Return the co-occurrence table underlying [`nestingplot`](@ref), in long
format: one row per combination of a level of grouping factor `A` and a level
of grouping factor `B`, for every pair of distinct grouping factors in `m`
(restricted/ordered by `gfs`, as in `nestingplot`). Columns:
- `factor_a`, `level_a`, `factor_b`, `level_b`: the pair of levels
- `count`: the number of observations sharing that pair of levels (`0` for
  combinations that never co-occur — these absences are what reveal nesting
  or partial crossing)

See [`nestingstructure`](@ref) for the pairwise nested/crossed classification
instead of the raw counts.
"""
function nestingtable(m::MixedModel; gfs::Union{Nothing,AbstractVector}=nothing)
    return _nesting_incidence_table(_nesting_data(m; gfs))
end

"""
    _nesting_structure_table(info::NamedTuple)

Convert the `info` NamedTuple (as returned by [`_nesting_data`](@ref)) into a
`DataFrame`: one row per pair of distinct grouping factors, with the pairwise
classification computed by [`_nesting_relationship`](@ref).
"""
function _nesting_structure_table(info::NamedTuple)
    factor_a = String[]
    factor_b = String[]
    n_levels_a = Int[]
    n_levels_b = Int[]
    relationship = Symbol[]
    density = Float64[]

    for (i, j) in sort!(collect(keys(info.relationships)))
        rel = info.relationships[(i, j)]
        push!(factor_a, info.names[i])
        push!(factor_b, info.names[j])
        push!(n_levels_a, info.nlevels[i])
        push!(n_levels_b, info.nlevels[j])
        push!(relationship, rel.rel)
        push!(density, rel.density)
    end

    return DataFrame(; factor_a, factor_b, n_levels_a, n_levels_b, relationship,
                     density)
end

"""
    nestingstructure(m::MixedModel; gfs::Union{Nothing,AbstractVector}=nothing)

Return the pairwise nested/crossed classification underlying [`nestingplot`](@ref)'s
upper-triangle badges: one row per pair of distinct grouping factors in `m`
(restricted/ordered by `gfs`, as in `nestingplot`), with:
- `factor_a`, `factor_b`: the pair of grouping factors
- `n_levels_a`, `n_levels_b`: their numbers of levels
- `relationship`: one of `:identical`, `:a_nested_in_b`, `:b_nested_in_a`,
  `:complete_crossing`, or `:partial_crossing` (see
  [`_nesting_relationship`](@ref) for the exact definitions)
- `density`: the fraction of the `n_levels_a × n_levels_b` grid of
  combinations that is actually observed (`1.0` for `:complete_crossing`)

See [`nestingtable`](@ref) for the raw per-level co-occurrence counts instead
of this per-pair summary.
"""
function nestingstructure(m::MixedModel; gfs::Union{Nothing,AbstractVector}=nothing)
    return _nesting_structure_table(_nesting_data(m; gfs))
end

"""
    _relationship_label(rel::NamedTuple, name_a::AbstractString, name_b::AbstractString)

Short text label describing the relationship classification `rel` (as
returned by [`_nesting_relationship`](@ref)) between row-factor `name_a` and
column-factor `name_b`.
"""
function _relationship_label(rel::NamedTuple, name_a::AbstractString,
                             name_b::AbstractString)
    txt = if rel.rel === :identical
        "identical"
    elseif rel.rel === :b_nested_in_a
        "$(name_b) ⊂ $(name_a)"
    elseif rel.rel === :a_nested_in_b
        "$(name_a) ⊂ $(name_b)"
    elseif rel.rel === :complete_crossing
        "crossed\n(complete)"
    else
        "crossed\n(partial)"
    end
    if rel.rel === :partial_crossing
        txt *= "\n$(round(100 * rel.density; digits=1))%"
    end
    return txt
end

"""
    _nestingplot_render!(f::Indexable, info::NamedTuple; colormap=:Blues, fontsize=14)

Render the grouping-factor nesting/crossing matrix into `f` from a
pre-computed `info` NamedTuple (as returned by [`_nesting_data`](@ref)).
"""
function _nestingplot_render!(f::Indexable, info::NamedTuple;
                              colormap=:Blues, fontsize::Real=14)
    k = length(info.names)
    gl = GridLayout()
    f[1, 1] = gl
    for i in 1:k, j in 1:k
        ax = Axis(gl[i, j])
        if i == j
            text!(ax, 0.5, 0.6; text=info.names[i], align=(:center, :center),
                  space=:relative, fontsize)
            text!(ax, 0.5, 0.35; text="n = $(info.nlevels[i])",
                  align=(:center, :center), space=:relative, fontsize=fontsize - 2,
                  color=:gray)
            limits!(ax, 0, 1, 0, 1)
            hidedecorations!(ax)
        elseif i > j
            tab = info.tables[(i, j)]
            ro, co = _block_order(tab)
            heatmap!(ax, permutedims(tab[ro, co]); colormap)
            hidedecorations!(ax; grid=false)
        else # i < j: mirror of the (j, i) pair, drawn as a text badge
            rel = info.relationships[(j, i)]
            label = _relationship_label(rel, info.names[j], info.names[i])
            text!(ax, 0.5, 0.5; text=label, align=(:center, :center),
                  space=:relative, fontsize=fontsize - 2)
            limits!(ax, 0, 1, 0, 1)
            hidedecorations!(ax)
        end
    end
    for i in 1:k
        rowsize!(gl, i, Relative(1 / k))
        colsize!(gl, i, Relative(1 / k))
    end
    colgap!(gl, 0)
    rowgap!(gl, 0)
    return f
end

"""
    nestingplot!(f::Indexable, m::MixedModel;
                 gfs::Union{Nothing,AbstractVector}=nothing,
                 colormap=:Blues, fontsize::Real=14)

Add a nesting/crossing matrix for the grouping factors of `m` to `f`.

Each distinct grouping factor (e.g. `subj`, `item`) becomes one row/column of
a matrix, laid out like a correlation-matrix plot:
- **diagonal**: factor name and number of levels
- **lower triangle**: a heatmap of the co-occurrence contingency table between
  the row and column factor. Rows and columns are reordered (for display only)
  to group each level with its most common partner, so that nested structure
  appears as a block-diagonal pattern and crossed structure appears as a dense
  or scattered rectangle.
- **upper triangle**: a text badge classifying the same pair as one factor
  nested in the other (`A ⊂ B`), `identical` (the two names partition
  observations the same way), or `crossed` — `complete` if every combination
  of levels is observed, `partial` (with the observed density) otherwise.

`gfs`, if given, restricts and orders which grouping factors to show (by name
or index); otherwise all distinct grouping factors are shown, in the order
they first appear in the model.

The classification only depends on the grouping factors' level assignments
(`ReMat.refs`), not on the original data, matching the model-based approach
used by [`upsetplot`](@ref).
"""
function nestingplot!(f::Indexable, m::MixedModel;
                      gfs::Union{Nothing,AbstractVector}=nothing,
                      colormap=:Blues, fontsize::Real=14)
    info = _nesting_data(m; gfs)
    return _nestingplot_render!(f, info; colormap, fontsize)
end

"""
    nestingplot(m::MixedModel;
                gfs::Union{Nothing,AbstractVector}=nothing,
                colormap=:Blues, fontsize::Real=14)

Return a `Figure` with a nesting/crossing matrix for the grouping factors of
`m`. See [`nestingplot!`](@ref) for details.
"""
function nestingplot(m::MixedModel; kwargs...)
    return nestingplot!(Figure(; size=(800, 800)), m; kwargs...)
end
