"""
    _categorical_terms(m::MixedModel) → (names, levels)

Return parallel vectors of categorical fixed-effect predictor names and their
level vectors (all levels, including the reference level).

Uses formula introspection: `first(m.formula.rhs)` is the `MatrixTerm` for
fixed effects; `StatsModels.terms` recursively collects leaf terms.
"""
function _categorical_terms(m::MixedModel)
    fe = first(m.formula.rhs)
    leaf = terms(fe)
    names = String[]
    lvls = Vector[]
    for term in leaf
        term isa CategoricalTerm || continue
        n = string(term.sym)
        n in names && continue
        push!(names, n)
        push!(lvls, term.contrasts.levels)
    end
    return names, lvls
end

"""
    _obs_predictor_cols(m::MixedModel) → Vector{Vector{String}}

Recover per-observation categorical predictor values from the model's fixed-effects
design matrix and contrast coding. Returns one `Vector{String}` per predictor,
in the same order as `_categorical_terms`.
"""
function _obs_predictor_cols(m::MixedModel)
    cn = coefnames(m)
    fe = first(m.formula.rhs)
    seen = Set{Symbol}()
    pred_cols = Vector{String}[]

    for term in terms(fe)
        term isa CategoricalTerm || continue
        term.sym in seen && continue
        push!(seen, term.sym)

        cm = term.contrasts.matrix
        lvls = term.contrasts.levels
        tnames = [string(term.sym, ": ", tn) for tn in term.contrasts.termnames]
        col_idxs = [findfirst(==(tn), cn) for tn in tnames]

        n_obs = size(m.X, 1)
        col = Vector{String}(undef, n_obs)
        for obs_i in 1:n_obs
            for (li, lv) in enumerate(lvls)
                if view(m.X, obs_i, col_idxs) ≈ view(cm, li, :)
                    col[obs_i] = string(lv)
                    break
                end
            end
        end
        push!(pred_cols, col)
    end

    return pred_cols
end

"""
    _upset_core(pred_names, pred_levels, pred_cols, gf_refs, gf_levels)

Shared computation for UpSet plots: build sets, full factorial cells, structural
combination matrix, counts, and marginal cells from already-identified predictor
names/levels. Called by both `_upset_data` (model path) and
`_upset_data_from_table` (table path).

`pred_cols` is a vector of string vectors, one per predictor, giving the
per-observation level. `gf_refs` is either `nothing` (count observations) or
an integer vector mapping each observation to a grouping-factor level index.

Returns a `NamedTuple` with fields:
`gf_levels`, `set_labels`, `cell_labels`, `combo_matrix`, `cell_counts`,
`cell_degrees`, `set_counts`.
"""
function _upset_core(pred_names::Vector{String}, pred_levels::Vector,
                     pred_cols::Vector{Vector{String}},
                     gf_refs::Union{Vector{<:Integer},Nothing},
                     gf_levels::Union{Vector,Nothing})
    # Sets = individual condition levels, grouped by predictor (columns of matrix)
    set_labels = String[string(p, ": ", lv)
                        for (p, lvs) in zip(pred_names, pred_levels) for lv in lvs]
    n_sets = length(set_labels)

    # Index: (predictor_index, level_as_string) → set column index
    set_index = Dict{Tuple{Int,String},Int}()
    k = 1
    for (pi, lvs) in enumerate(pred_levels)
        for lv in lvs
            set_index[(pi, string(lv))] = k
            k += 1
        end
    end

    n_expected = prod(length.(pred_levels))
    if n_expected > 1000
        dims = join(["$(p) ($(length(l)))" for (p, l) in zip(pred_names, pred_levels)],
                    " × ")
        throw(ArgumentError("Full factorial has $n_expected cells ($dims); the plot would be unreadable. " *
                            "Use `cols` to select fewer or lower-cardinality predictors."))
    end

    cell_labels = String[]
    cell_set_indices = Vector{Int}[]
    cell_combo_strs = Vector{String}[]

    for combo in Iterators.product(pred_levels...)
        strs = [string(combo[j]) for j in eachindex(pred_names)]
        parts = [string(pred_names[j], ": ", strs[j]) for j in eachindex(pred_names)]
        push!(cell_labels, join(parts, " & "))
        push!(cell_set_indices, [set_index[(j, strs[j])] for j in eachindex(pred_names)])
        push!(cell_combo_strs, strs)
    end
    n_cells = length(cell_labels)
    cell_lookup = Dict(Tuple(strs) => ci for (ci, strs) in enumerate(cell_combo_strs))

    # Structural combination matrix: (n_cells × n_sets)
    combo_matrix = falses(n_cells, n_sets)
    for (ci, sidxs) in enumerate(cell_set_indices)
        for si in sidxs
            combo_matrix[ci, si] = true
        end
    end

    n_obs = length(first(pred_cols))

    if gf_refs !== nothing
        n_gf = length(gf_levels)
        cell_membership = falses(n_gf, n_cells)
        set_membership = falses(n_gf, n_sets)

        for obs_i in 1:n_obs
            gi = Int(gf_refs[obs_i])
            obs_vals = Tuple(pred_cols[j][obs_i] for j in eachindex(pred_names))

            ci = get(cell_lookup, obs_vals, nothing)
            if ci !== nothing
                cell_membership[gi, ci] = true
            end

            for j in eachindex(pred_names)
                si = get(set_index, (j, obs_vals[j]), nothing)
                isnothing(si) && continue
                set_membership[gi, si] = true
            end
        end

        cell_counts = [count(cell_membership[:, ci]) for ci in axes(cell_membership, 2)]
        set_counts = [count(set_membership[:, si]) for si in axes(set_membership, 2)]
    else
        cell_counts = zeros(Int, n_cells)
        set_counts = zeros(Int, n_sets)

        for obs_i in 1:n_obs
            obs_vals = Tuple(pred_cols[j][obs_i] for j in eachindex(pred_names))

            ci = get(cell_lookup, obs_vals, nothing)
            if ci !== nothing
                cell_counts[ci] += 1
            end

            for j in eachindex(pred_names)
                si = get(set_index, (j, obs_vals[j]), nothing)
                isnothing(si) && continue
                set_counts[si] += 1
            end
        end
    end

    n_preds = length(pred_names)
    cell_degrees = fill(n_preds, n_cells)

    # Marginal cells: collapse exactly one predictor (all its levels active).
    marginal_labels = String[]
    marginal_combo_rows = Vector{Int}[]
    marginal_counts = Int[]
    marginal_degrees = Int[]

    for collapse_pi in eachindex(pred_names)
        other_pis = [j for j in eachindex(pred_names) if j != collapse_pi]
        other_level_seqs = [pred_levels[j] for j in other_pis]

        for other_combo in Iterators.product(other_level_seqs...)
            other_strs = [string(other_combo[k]) for k in eachindex(other_pis)]

            group = [ci
                     for (ci, cv) in enumerate(cell_combo_strs)
                     if all(cv[other_pis[k]] == other_strs[k] for k in eachindex(other_pis))]

            label_parts = [string(pred_names[other_pis[k]], ": ", other_strs[k])
                           for k in eachindex(other_pis)]
            push!(marginal_labels, join(label_parts, " & "))

            filled = Int[]
            for lv in pred_levels[collapse_pi]
                push!(filled, set_index[(collapse_pi, string(lv))])
            end
            for k in eachindex(other_pis)
                push!(filled, set_index[(other_pis[k], other_strs[k])])
            end
            push!(marginal_combo_rows, sort!(filled))

            if gf_refs !== nothing
                m_mem = trues(length(gf_levels))
                for ci in group
                    m_mem .&= cell_membership[:, ci]
                end
                push!(marginal_counts, count(m_mem))
            else
                push!(marginal_counts, sum(cell_counts[ci] for ci in group; init=0))
            end

            push!(marginal_degrees, n_preds - 1)
        end
    end

    n_marginals = length(marginal_labels)
    marginal_combo = falses(n_marginals, n_sets)
    for (mi, filled) in enumerate(marginal_combo_rows)
        for si in filled
            marginal_combo[mi, si] = true
        end
    end

    all_labels = vcat(cell_labels, marginal_labels)
    all_combo = vcat(combo_matrix, marginal_combo)
    all_counts = vcat(cell_counts, marginal_counts)
    all_degrees = vcat(cell_degrees, marginal_degrees)

    return (; gf_levels, set_labels, cell_labels=all_labels, combo_matrix=all_combo,
            cell_counts=all_counts, cell_degrees=all_degrees, set_counts)
end

"""
    _upset_data(m::MixedModel, gf=first(fnames(m)))

Build UpSet data structures from a fitted model.

Predictor names, levels, and per-observation values are recovered from the
model's formula, design matrix, and contrast coding.
Grouping-factor membership comes from the random-effects terms.
"""
function _upset_data(m::MixedModel, gf::Union{Symbol,Nothing}=first(fnames(m)))
    pred_names, pred_levels = _categorical_terms(m)
    isempty(pred_names) &&
        throw(ArgumentError("No categorical fixed-effect predictors found in model"))

    pred_cols = _obs_predictor_cols(m)

    if gf !== nothing
        idx = findfirst(==(gf), fnames(m))
        isnothing(idx) &&
            throw(ArgumentError("$gf is not the name of a grouping variable in the model"))
        gf_refs = m.reterms[idx].refs
        gf_levels = m.reterms[idx].levels
    else
        gf_refs = nothing
        gf_levels = nothing
    end

    return _upset_core(pred_names, pred_levels, pred_cols, gf_refs, gf_levels)
end

"""
    _upset_data_from_table(data; cols=All(), gf=nothing)

Build UpSet data structures directly from a Tables.jl-compatible table.

Non-numeric columns (after applying the `cols` selector) are treated as
categorical predictors. The `gf` column, if specified, is excluded from the
predictor set and used as a grouping factor for counting unique levels.
"""
function _upset_data_from_table(data; cols=All(),
                                gf::Union{AbstractString,Symbol,Nothing}=nothing)
    df = DataFrame(data)
    gf = gf === nothing ? nothing : Symbol(gf)
    candidate_names = propertynames(select(df, cols))
    cat_names = [n
                 for n in candidate_names
                 if !(nonmissingtype(eltype(df[!, n])) <: Number) && n !== gf]
    isempty(cat_names) &&
        throw(ArgumentError("No categorical columns found in data after filtering"))

    pred_names = string.(cat_names)
    pred_levels = [sort!(unique(string.(skipmissing(df[!, n])))) for n in cat_names]

    for (n, lvs) in zip(cat_names, pred_levels)
        if length(lvs) > 10
            @warn "Column $n has $(length(lvs)) distinct levels and may be an " *
                  "identifier column rather than a categorical predictor. " *
                  "Consider excluding it via `cols`."
        end
    end

    pred_cols = [string.(df[!, n]) for n in cat_names]

    if gf !== nothing
        gf_levels = sort!(unique(df[!, gf]))
        gf_index = Dict(lv => i for (i, lv) in enumerate(gf_levels))
        gf_refs = [gf_index[v] for v in df[!, gf]]
    else
        gf_refs = nothing
        gf_levels = nothing
    end

    return _upset_core(pred_names, pred_levels, pred_cols, gf_refs, gf_levels)
end

"""
    _upset_incidence_table(info::NamedTuple)

Convert the `info` NamedTuple (as returned by [`_upset_data`](@ref) or
[`_upset_data_from_table`](@ref)) into a `DataFrame`: one row per combination
cell, with `cell`, `degree`, and `count` columns, plus one `Bool` column per
set giving its raw label (e.g. `"gender: M"`).
"""
function _upset_incidence_table(info::NamedTuple)
    df = DataFrame(; cell=info.cell_labels, degree=info.cell_degrees,
                   count=info.cell_counts)
    for (si, label) in enumerate(info.set_labels)
        df[!, Symbol(label)] = info.combo_matrix[:, si]
    end
    return df
end

"""
    upsettable(m::MixedModel, gf::Union{Symbol,Nothing}=first(fnames(m)))

Return the incidence table underlying [`upsetplot`](@ref): one row per
combination cell (every full factorial cell of the categorical fixed-effect
predictors, plus marginal cells that collapse a single predictor), with:
- `cell`: a label for the combination
- `degree`: the number of predictors involved (one fewer for marginal cells)
- `count`: the number of observations (or grouping-factor levels, if `gf` is
  given) falling in the cell
- one `Bool` column per set (an individual condition level, e.g.
  `"gender: M"`) indicating whether that set is active in the cell

Pass `gf=nothing` to count observations instead of grouping-factor levels,
matching [`upsetplot`](@ref).
"""
function upsettable(m::MixedModel, gf::Union{Symbol,Nothing}=first(fnames(m)))
    return _upset_incidence_table(_upset_data(m, gf))
end

"""
    upsettable(data; cols=All(), gf::Union{Symbol,Nothing}=nothing)

Return the incidence table underlying [`upsetplot`](@ref), built directly
from a Tables.jl-compatible table. See [`upsetplot`](@ref) for the meaning of
`cols` and `gf`.
"""
function upsettable(data; cols=All(), gf::Union{Symbol,Nothing}=nothing)
    return _upset_incidence_table(_upset_data_from_table(data; cols, gf))
end

"""
    _upsetplot_render!(f::Indexable, info::NamedTuple; kwargs...)

Render the UpSet combination matrix, intersection bar chart, and (optionally)
set-size bar chart into `f` from a pre-computed `info` NamedTuple (as returned
by `_upset_data` or `_upset_data_from_table`).

Layout keywords:
- `orientation`: `:horizontal` (combinations as columns, sets as rows — default)
  or `:vertical` (combinations as rows, sets as columns).
- `show_setsize`: whether to show the set-size bar chart (`true` by default).
- `intersection_pos`: position of the intersection-size bars relative to the
  incidence matrix — `:top`/`:bottom` for horizontal, `:left`/`:right` for
  vertical.  Defaults to `:top` (horizontal) or `:left` (vertical).
- `setsize_pos`: position of the set-size bars — `:left`/`:right` for
  horizontal, `:top`/`:bottom` for vertical.  Defaults to `:left` (horizontal)
  or `:top` (vertical).
"""
function _upsetplot_render!(f::Indexable, info::NamedTuple;
                            sortby::Symbol=:count,
                            show_empty::Bool=true,
                            orientation::Symbol=:horizontal,
                            show_setsize::Bool=true,
                            intersection_pos::Union{Symbol,Nothing}=nothing,
                            setsize_pos::Union{Symbol,Nothing}=nothing,
                            filled_color=:black,
                            empty_color=:lightgray,
                            bar_color=:steelblue,
                            dot_size=12)
    orientation in (:horizontal, :vertical) ||
        throw(ArgumentError("orientation must be :horizontal or :vertical, got :$orientation"))
    horiz = orientation === :horizontal

    if intersection_pos === nothing
        intersection_pos = horiz ? :top : :left
    end
    if setsize_pos === nothing
        setsize_pos = horiz ? :left : :top
    end

    if horiz
        intersection_pos in (:top, :bottom) ||
            throw(ArgumentError("intersection_pos must be :top or :bottom for horizontal orientation, got :$intersection_pos"))
        if show_setsize
            setsize_pos in (:left, :right) ||
                throw(ArgumentError("setsize_pos must be :left or :right for horizontal orientation, got :$setsize_pos"))
        end
    else
        intersection_pos in (:left, :right) ||
            throw(ArgumentError("intersection_pos must be :left or :right for vertical orientation, got :$intersection_pos"))
        if show_setsize
            setsize_pos in (:top, :bottom) ||
                throw(ArgumentError("setsize_pos must be :top or :bottom for vertical orientation, got :$setsize_pos"))
        end
    end

    n_sets = length(info.set_labels)

    perm = if sortby === :count
        sortperm(info.cell_counts; rev=true)
    elseif sortby === :degree
        sortperm(collect(zip(info.cell_degrees, info.cell_counts)); rev=true)
    else
        throw(ArgumentError("sortby must be :count or :degree, got :$sortby"))
    end
    perm = show_empty ? perm : filter(i -> info.cell_counts[i] > 0, perm)
    cell_counts = info.cell_counts[perm]
    combo_matrix = info.combo_matrix[perm, :]
    n_shown = length(cell_counts)

    # --- Axes: grid positions, bar charts, decorations, linking ---
    if horiz
        bar_row = intersection_pos === :top ? 1 : 2
        matrix_row = intersection_pos === :top ? 2 : 1
        if show_setsize
            sets_col = setsize_pos === :left ? 1 : 2
            matrix_col = setsize_pos === :left ? 2 : 1
        else
            matrix_col = 1
        end

        ax_bar = Axis(f[bar_row, matrix_col]; ylabel="Intersection size")
        ax_matrix = Axis(f[matrix_row, matrix_col])
        hidexdecorations!(ax_bar; grid=false)
        hidexdecorations!(ax_matrix; grid=false)
        linkxaxes!(ax_bar, ax_matrix)

        barplot!(ax_bar, 1:n_shown, cell_counts; color=bar_color)
        if intersection_pos === :bottom
            ax_bar.yreversed = true
        end

        ax_matrix.yticks = (1:n_sets, info.set_labels)
        ax_matrix.yreversed = true

        if show_setsize
            ax_sets = Axis(f[matrix_row, sets_col]; xlabel="Set size")
            hideydecorations!(ax_sets; grid=false)
            linkyaxes!(ax_sets, ax_matrix)
            barplot!(ax_sets, 1:n_sets, info.set_counts; direction=:x, color=bar_color)
            ax_sets.yreversed = true
            if setsize_pos === :left
                ax_sets.xreversed = true
            end
        end
    else
        bar_col = intersection_pos === :left ? 1 : 2
        matrix_col = intersection_pos === :left ? 2 : 1
        if show_setsize
            sets_row = setsize_pos === :top ? 1 : 2
            matrix_row = setsize_pos === :top ? 2 : 1
        else
            matrix_row = 1
        end

        ax_bar = Axis(f[matrix_row, bar_col]; xlabel="Intersection size")
        ax_matrix = Axis(f[matrix_row, matrix_col])
        hideydecorations!(ax_bar; grid=false)
        hideydecorations!(ax_matrix; grid=false)
        linkyaxes!(ax_bar, ax_matrix)

        barplot!(ax_bar, 1:n_shown, cell_counts; direction=:x, color=bar_color)
        ax_bar.yreversed = true
        if intersection_pos === :left
            ax_bar.xreversed = true
        end

        ax_matrix.xticks = (1:n_sets, info.set_labels)
        ax_matrix.xticklabelrotation = π / 4
        ax_matrix.yreversed = true

        if show_setsize
            if setsize_pos === :bottom
                ax_matrix.xaxisposition = :top
            end
            ax_sets = Axis(f[sets_row, matrix_col]; ylabel="Set size")
            hidexdecorations!(ax_sets; grid=false)
            linkxaxes!(ax_sets, ax_matrix)
            barplot!(ax_sets, 1:n_sets, info.set_counts; color=bar_color)
            if setsize_pos === :bottom
                ax_sets.yreversed = true
            end
        end
    end

    # --- Combination matrix: dots and connecting lines ---
    empty_combo = Float64[]
    empty_set = Float64[]
    filled_combo = Float64[]
    filled_set = Float64[]

    for ci in 1:n_shown
        active = findall(combo_matrix[ci, :])
        if length(active) >= 2
            s_min, s_max = extrema(active)
            if horiz
                lines!(ax_matrix, [ci, ci], [s_min, s_max];
                       color=filled_color, linewidth=2)
            else
                lines!(ax_matrix, [s_min, s_max], [ci, ci];
                       color=filled_color, linewidth=2)
            end
        end
        for si in 1:n_sets
            if combo_matrix[ci, si]
                push!(filled_combo, ci)
                push!(filled_set, si)
            else
                push!(empty_combo, ci)
                push!(empty_set, si)
            end
        end
    end

    if horiz
        ex, ey = empty_combo, empty_set
        fx, fy = filled_combo, filled_set
    else
        ex, ey = empty_set, empty_combo
        fx, fy = filled_set, filled_combo
    end

    isempty(ex) ||
        scatter!(ax_matrix, ex, ey; color=empty_color, markersize=dot_size)
    isempty(fx) ||
        scatter!(ax_matrix, fx, fy; color=filled_color, markersize=dot_size)

    return f
end

"""
    upsetplot!(f::Indexable, m::MixedModel,
               gf::Union{Symbol,Nothing}=first(fnames(m));
               sortby::Symbol=:count,
               show_empty::Bool=true,
               orientation::Symbol=:horizontal,
               show_setsize::Bool=true,
               intersection_pos::Symbol=...,
               setsize_pos::Symbol=...,
               filled_color=:black,
               empty_color=:lightgray,
               bar_color=:steelblue,
               dot_size=12)

Add an UpSet plot to `f` showing which levels of grouping factor `gf` appear
in which categorical fixed-effect conditions.

Predictor names, levels, and per-observation values are recovered from the
model's formula, design matrix, and contrast coding — no original data frame
is needed.

**Layout keywords:**
- `orientation`: `:horizontal` (default — combinations as columns, sets as rows)
  or `:vertical` (combinations as rows, sets as columns).
- `show_setsize`: show the set-size bar chart (`true` by default).
- `intersection_pos`: where the intersection-size bars go relative to the
  incidence matrix — `:top`/`:bottom` (horizontal) or `:left`/`:right`
  (vertical).
- `setsize_pos`: where the set-size bars go — `:left`/`:right` (horizontal) or
  `:top`/`:bottom` (vertical).
"""
function upsetplot!(f::Indexable, m::MixedModel,
                    gf::Union{Symbol,Nothing}=first(fnames(m));
                    kwargs...)
    info = _upset_data(m, gf)
    return _upsetplot_render!(f, info; kwargs...)
end

"""
    upsetplot!(f::Indexable, data;
               cols=All(),
               gf::Union{Symbol,Nothing}=nothing,
               sortby::Symbol=:count,
               show_empty::Bool=true,
               orientation::Symbol=:horizontal,
               show_setsize::Bool=true,
               intersection_pos::Symbol=...,
               setsize_pos::Symbol=...,
               filled_color=:black,
               empty_color=:lightgray,
               bar_color=:steelblue,
               dot_size=12)

Add an UpSet plot to `f` directly from a Tables.jl-compatible table.

Non-numeric columns (optionally restricted by `cols`) become the sets. Pass
`gf=:col` to count unique values of that column per cell instead of
observations.

See [`upsetplot!(::Indexable, ::MixedModel)`](@ref) for layout keyword details.
"""
function upsetplot!(f::Indexable, data;
                    cols=All(),
                    gf::Union{Symbol,Nothing}=nothing,
                    kwargs...)
    info = _upset_data_from_table(data; cols, gf)
    return _upsetplot_render!(f, info; kwargs...)
end

"""
    upsetplot(m::MixedModel,
              gf::Union{Symbol,Nothing}=first(fnames(m));
              kwargs...)

Return a `Figure` with an UpSet plot showing which levels of grouping factor `gf`
appear in which categorical fixed-effect conditions.

Pass `gf=nothing` to count observations instead of grouping-factor levels.

Predictor names, levels, and per-observation values are recovered from the
model's formula, design matrix, and contrast coding — no original data frame
is needed.

`kwargs` are forwarded to [`upsetplot!`](@ref).
"""
function upsetplot(m::MixedModel, args...; kwargs...)
    return upsetplot!(Figure(; size=(1000, 800)), m, args...; kwargs...)
end

"""
    upsetplot(data; cols=All(), gf::Union{Symbol,Nothing}=nothing, kwargs...)

Return a `Figure` with an UpSet plot built directly from a Tables.jl-compatible
table.

Non-numeric columns (optionally restricted by `cols`) are used as sets. Pass
`gf=:col` to count unique values of that column per cell instead of observations.

`kwargs` are forwarded to [`upsetplot!`](@ref).
"""
function upsetplot(data; cols=All(), gf::Union{Symbol,Nothing}=nothing, kwargs...)
    return upsetplot!(Figure(; size=(1000, 800)), data; cols, gf, kwargs...)
end
