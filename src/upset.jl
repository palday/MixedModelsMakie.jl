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
function _upset_core(pred_names, pred_levels, pred_cols::AbstractVector,
                     gf_refs::Union{AbstractVector{<:Integer},Nothing}, gf_levels)
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
model's formula, design matrix, and contrast coding — no original data frame
is needed. Grouping-factor membership comes from the random-effects terms.
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
    _upsetplot_render!(f::Indexable, info::NamedTuple; kwargs...)

Render the UpSet combination matrix, intersection bar chart, and set-size bar
chart into `f` from a pre-computed `info` NamedTuple (as returned by
`_upset_data` or `_upset_data_from_table`).
"""
function _upsetplot_render!(f::Indexable, info::NamedTuple;
                            sortby::Symbol=:count,
                            show_empty::Bool=true,
                            filled_color=:black,
                            empty_color=:lightgray,
                            bar_color=:steelblue,
                            dot_size=12)
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

    ax_bar = Axis(f[1, 2]; ylabel="Intersection size")
    ax_matrix = Axis(f[2, 2])
    ax_sets = Axis(f[2, 1]; xlabel="Set size")

    hidexdecorations!(ax_bar; grid=false)
    hidexdecorations!(ax_matrix; grid=false)
    hideydecorations!(ax_sets; grid=false)
    linkxaxes!(ax_bar, ax_matrix)
    linkyaxes!(ax_sets, ax_matrix)

    barplot!(ax_bar, 1:n_shown, cell_counts; color=bar_color)

    barplot!(ax_sets, 1:n_sets, info.set_counts; direction=:x, color=bar_color)
    ax_sets.xreversed = true

    empty_xs = Float64[]
    empty_ys = Float64[]
    filled_xs = Float64[]
    filled_ys = Float64[]

    for ci in 1:n_shown
        active = findall(combo_matrix[ci, :])
        if length(active) >= 2
            lines!(ax_matrix, [ci, ci], [minimum(active), maximum(active)];
                   color=filled_color, linewidth=2)
        end
        for si in 1:n_sets
            if combo_matrix[ci, si]
                push!(filled_xs, ci)
                push!(filled_ys, si)
            else
                push!(empty_xs, ci)
                push!(empty_ys, si)
            end
        end
    end

    isempty(empty_xs) ||
        scatter!(ax_matrix, empty_xs, empty_ys; color=empty_color, markersize=dot_size)
    isempty(filled_xs) ||
        scatter!(ax_matrix, filled_xs, filled_ys; color=filled_color, markersize=dot_size)

    ax_matrix.yticks = (1:n_sets, info.set_labels)
    ax_matrix.yreversed = true
    ax_sets.yreversed = true

    return f
end

"""
    upsetplot!(f::Indexable, m::MixedModel,
               gf::Union{Symbol,Nothing}=first(fnames(m));
               sortby::Symbol=:count,
               show_empty::Bool=true,
               filled_color=:black,
               empty_color=:lightgray,
               bar_color=:steelblue,
               dot_size=12)

Add an UpSet plot to `f` showing which levels of grouping factor `gf` appear
in which categorical fixed-effect conditions.

Predictor names, levels, and per-observation values are recovered from the
model's formula, design matrix, and contrast coding — no original data frame
is needed.

**Layout (standard UpSet orientation):**
- top-right: intersection-size bar chart (e.g., subjects per full factorial cell)
- middle-right: combination matrix — rows are individual condition levels (sets),
  columns are full factorial cells; filled circles mark which condition levels are
  active in each cell, connected by a vertical line
- middle-left: set-size bar chart (e.g., subjects per individual condition level)
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
               filled_color=:black,
               empty_color=:lightgray,
               bar_color=:steelblue,
               dot_size=12)

Add an UpSet plot to `f` directly from a Tables.jl-compatible table.

Non-numeric columns (optionally restricted by `cols`) become the sets. Pass
`gf=:col` to count unique values of that column per cell instead of
observations.
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
