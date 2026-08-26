function _coefnames(x::MixedModel, ptype::Nothing=nothing; show_intercept=true,
                    group=nothing)
    isnothing(group) || throw(ArgumentError("`group` not supported for MixedModel"))
    cn = fixefnames(x)
    return show_intercept ? cn : filter!(!=("(Intercept)"), cn)
end

"""
    _normalize_ptype(ptype)

Map the ASCII aliases `:sigma`, `:rho`, `:theta`, `:beta` to their Greek `ptype`
symbols (`:σ`, `:ρ`, `:θ`, `:β`); any other value (including `nothing` or `:β`)
passes through unchanged.
"""
function _normalize_ptype(ptype)
    ptype === :sigma && return :σ
    ptype === :rho && return :ρ
    ptype === :theta && return :θ
    ptype === :beta && return :β
    return ptype
end

"""
    _validate_group(ptype, group)

Throw an `ArgumentError` if `group` is specified for a `ptype` that has no
associated grouping factor. Only `:σ`/`:ρ` are qualified by a grouping
factor in [`_allparlabel`](@ref)'s labeling scheme (`:β` has none, `:θ`'s
labels are purely positional).
"""
function _validate_group(ptype, group)
    group === nothing || ptype in (:σ, :ρ) ||
        throw(ArgumentError("`group` is only supported for ptype ∈ (:σ, :ρ); got ptype=$(ptype)"))
    return nothing
end

"""
    _parambounds(ptype, bsamp, coefname)

Return the `(lower, upper)` support bounds of a single bootstrap parameter,
used to truncate ridge density curves at their theoretically valid range
(kernel smoothing can otherwise leak density past a hard boundary, e.g.
negative σ or |ρ| > 1). `:β` is unbounded. `:σ` is bounded to `[0, Inf]`.
`:ρ` is bounded to `[-1, 1]`. `:θ` is bounded per-element via
`lowerbd(bsamp)` (upper bound `Inf`), matched to `coefname`'s positional
index (`θ01` → `lowerbd(bsamp)[1]`, etc.).
"""
_parambounds(ptype, bsamp, coefname) = _parambounds(Val(ptype), bsamp, coefname)
_parambounds(::Val{:β}, bsamp, coefname) = (-Inf, Inf)
_parambounds(::Val{:σ}, bsamp, coefname) = (0.0, Inf)
_parambounds(::Val{:ρ}, bsamp, coefname) = (-1.0, 1.0)
function _parambounds(::Val{:θ}, bsamp, coefname)
    idx = parse(Int, replace(coefname, "θ" => ""))
    return (lowerbd(bsamp)[idx], Inf)
end

"""
    _StepCurve(x, density)

A step-function outline with fields `x`/`density`, mirroring the shape of
`KernelDensity.UnivariateKDE` so [`_histcurve`](@ref)'s result can be drawn
with the same ridge-plotting code as a KDE. A plain struct (rather than a
`NamedTuple`) so `DataFrames.combine` stores it as a single cell instead of
trying to expand it into multiple columns.
"""
struct _StepCurve
    x::Vector{Float64}
    density::Vector{Float64}
end

"""
    _histcurve(vals; bins=nothing, bounds=(-Inf, Inf))

Return a [`_StepCurve`](@ref) outline of a histogram of `vals`. Unlike a
KDE, a histogram never draws mass *within* a bin past a parameter's hard
bounds. But `StatsBase`'s automatic bin edges are rounded to "nice" numbers
and can still overshoot the true data range, so the outermost edges are
clamped to `bounds` (see [`_parambounds`](@ref)) the same way a KDE curve
is truncated. A histogram can otherwise show a genuine spike/impulse when
many bootstrap draws land on a boundary (e.g. a singular fit).

`bins` is forwarded to `StatsBase.fit(Histogram, vals; nbins=bins)` when
given; otherwise `StatsBase`'s automatic bin selection is used.
"""
function _histcurve(vals; bins=nothing, bounds=(-Inf, Inf))
    h = isnothing(bins) ? fit(Histogram, vals) : fit(Histogram, vals; nbins=bins)
    edges = collect(only(h.edges))
    edges[1] = max(edges[1], first(bounds))
    edges[end] = min(edges[end], last(bounds))
    counts = h.weights
    n = length(counts)
    x = Vector{Float64}(undef, 2n + 2)
    density = Vector{Float64}(undef, 2n + 2)
    x[1] = edges[1]
    density[1] = 0.0
    for i in 1:n
        x[2i] = edges[i]
        density[2i] = counts[i]
        x[2i + 1] = edges[i + 1]
        density[2i + 1] = counts[i]
    end
    x[end] = edges[n + 1]
    density[end] = 0.0
    return _StepCurve(x, density)
end

"""
    _allparlabel(group, names)

Combine the `group`/`names` columns of `MixedModels.allpars` into a single
coefname string: bare `names` for fixed effects (`group === missing`),
`"residual"` for the residual standard deviation (`names === missing`), and
`"group: names"` (e.g. `"subj: (Intercept)"`) otherwise.
"""
_allparlabel(group::Missing, names) = string(names)
_allparlabel(group, names::Missing) = group
_allparlabel(group, names) = string(group, ": ", names)

"""
    _bootstrap_longtable(bsamp::MixedModelBootstrap, ptype; group=nothing)

Return a long-format DataFrame with columns `:iter`, `:coefname`, `:value`
for the requested `ptype` (`:β`, `:σ`, `:ρ`, or `:θ`).

`:β`/`:σ`/`:ρ` are drawn from `bsamp.allpars`, with `coefname` built via
[`_allparlabel`](@ref). `:θ` has no named/tidied accessor in MixedModels.jl,
so it is drawn from `bsamp.tbl` and labeled positionally (e.g. `θ01`, `θ02`,
...), matching the convention already used by [`ridge2d`](@ref).

`group` restricts `:σ`/`:ρ` to a single grouping factor (e.g. `:subj`); see
[`_validate_group`](@ref).
"""
function _bootstrap_longtable(bsamp::MixedModelBootstrap, ptype; group=nothing)
    _validate_group(ptype, group)
    if ptype === :θ
        tbl = bsamp.tbl
        θcols = collect(filter(startswith("θ"), string.(propertynames(tbl))))
        df = DataFrame(Tables.columntable(tbl))
        df.iter = 1:nrow(df)
        return stack(df, θcols; variable_name=:coefname, value_name=:value)
    else
        df = DataFrame(bsamp.allpars)
        filter!(:type => ==(string(ptype)), df)
        if group !== nothing
            filter!(:group => ==(string(group)), df)
            isempty(df) &&
                throw(ArgumentError("No $(ptype) parameters found for group $(group)."))
        end
        df.coefname = _allparlabel.(df.group, df.names)
        return select(df, :iter, :coefname, :value)
    end
end

function _coefnames(x::MixedModelBootstrap, ptype; show_intercept=true, group=nothing)
    ptype = _normalize_ptype(something(ptype, :β))
    _validate_group(ptype, group)
    if ptype === :β
        nt = first(x.fits).β
        cn = [string(k) for (k, v) in pairs(nt) if !isequal(v, -0.0)]
        return show_intercept ? cn : filter!(!=("(Intercept)"), cn)
    else
        # relies on `unique` preserving first-occurrence order so that
        # y-axis ticks come out in a stable, sensible order (grouping
        # factors/RE terms in the order MixedModels.jl emits them)
        return unique(_bootstrap_longtable(x, ptype; group).coefname)
    end
end

"""
    confint_table(x::MixedModel, level=0.95; show_intercept=true)
    confint_table(x::MixedModelBootstrap, level=0.95; show_intercept=true)

Return a DataFrame of coefficient names, point estimates and confidence intervals.

`level` specifies the confidence level, e.g. `0.95` for 95% confidence intervals.

For `MixedModels`, the intervals are computed using the standard errors and the asymptotic
Wald approximation (e.g. est±1.96*se for 95% intervals). For `MixedModelBootstrap`, the intervals
are computed using `shortestcovint`, and the point estimate re-computed as the mean of the bootstrap values.

For `MixedModelBootstrap`, `ptype` selects which parameters to summarize:
`:β` (fixed effects, default), `:σ` (random-effect standard deviations,
including the residual, labeled `"group: names"`/`"residual"`), `:ρ`
(random-effect correlations, labeled `"group: name_i, name_j"`), or `:θ`
(the unconstrained Cholesky parameterization, labeled positionally as
`θ01`, `θ02`, ...). `ptype` is not supported for a plain `MixedModel`. The
ASCII aliases `:sigma`, `:rho`, `:theta` are also accepted.
`show_intercept` only filters the fixed-effect `"(Intercept)"` row; it is a
no-op for `:σ`/`:ρ`/`:θ`.

`group` restricts `ptype ∈ (:σ, :ρ)` to a single grouping factor, e.g.
`group=:subj`. It is not supported for `:β`/`:θ` (which have no associated
grouping factor) or for a plain `MixedModel`.

The returned table has the following columns:
- `coefname`: the names of the coefficients
- `estimate`: the point estimates
- `lower`: the lower edge of the confidence interval
- `upper`: the upper edge of the confidence interval

!!! note
    This function is internal and may be removed in a future release
    without being considered a breaking change.
"""
function confint_table(x::MixedModel, level=0.95; ptype=nothing, show_intercept=true,
                       group=nothing)
    isnothing(ptype) || throw(ArgumentError("`ptype` not supported for MixedModel"))
    isnothing(group) || throw(ArgumentError("`group` not supported for MixedModel"))
    # taking from the lower tail
    semultiple = zquantile((1 - level) / 2)
    se = stderror(x)
    est = coef(x)

    df = DataFrame(;
                   coefname=coefnames(x),
                   estimate=coef(x),
                   # signs are 'swapped' b/c semultiple comes from the lower tail
                   lower=est + semultiple * se,
                   upper=est - semultiple * se)
    return filter!(:coefname => in(_coefnames(x; show_intercept)), df)
end

function confint_table(x::MixedModelBootstrap, level=0.95; ptype=:β, show_intercept=true,
                       group=nothing)
    ptype = _normalize_ptype(something(ptype, :β))
    ptype in (:β, :σ, :ρ, :θ) || throw(ArgumentError("ptype $(ptype) not supported"))
    _validate_group(ptype, group)
    if ptype === :θ
        # allpars (and hence shortestcovint(bsamp, level)) has no θ rows,
        # so θ needs its own grouped HDI computation
        long = _bootstrap_longtable(x, :θ)
        df = combine(groupby(long, :coefname),
                     :value => mean => :estimate,
                     :value => (v -> NamedTuple{(:lower, :upper)}(shortestcovint(v, level))) => [:lower,
                                                                                                 :upper])
    else
        hdi_df = DataFrame(shortestcovint(x, level))
        filter!(:type => ==(string(ptype)), hdi_df)
        group === nothing || filter!(:group => ==(string(group)), hdi_df)
        hdi_df.coefname = _allparlabel.(hdi_df.group, hdi_df.names)
        est = combine(groupby(_bootstrap_longtable(x, ptype; group), :coefname),
                      :value => mean => :estimate)
        df = innerjoin(est, select(hdi_df, :coefname, :lower, :upper); on=:coefname)
    end
    return filter!(:coefname => in(_coefnames(x, ptype; show_intercept, group)), df)
end

_npreds(args...; kwargs...) = length(_coefnames(args...; kwargs...))

"""
    _extract_title!(ax::Axis, kwargs)::Base.Pairs

If a title is present in kwargs, use it to set the axis title.

Returns kwargs without an entry for `title`.
"""
function _extract_title!(ax::Axis, kwargs)::Base.Pairs
    if :title in keys(kwargs)
        ax.title = kwargs[:title]
        kwargs = NamedTuple((k => v for (k, v) in kwargs if k != :title))
    end
    return Base.pairs(kwargs)
end

function _place_legend!(figure, axis, position; kwargs...)
    if position === true
        position = :bottom
    elseif position === false
        return figure
    end
    if position === :bottom || position === :top
        orientation = :horizontal
        x = position === :top ? 0 : 2
        y = 1
    elseif position === :left || position === :right
        orientation = :vertical
        x = 1
        y = position === :left ? 0 : 2
    else
        throw(ArgumentError("Invalid legend position: $position"))
    end
    figure[x, y] = Legend(figure, axis;
                          orientation,
                          tell_width=false,
                          tell_height=false,
                          kwargs...)
    return figure
end

"""
    ppoints(n::Integer)

Return a sequence of `n` equally-spaced points in the interval (0, 1) - so-called "probability points"
"""
ppoints(n::Integer) = inv(2n):inv(n):1

"""
    zquantile(x::Real)

Evaluate `quantile(Normal(), x)` using only the `SpecialFunctions` package (i.e. not requiring `Distributions`).
"""
zquantile(x::T) where {T<:Real} = -erfcinv(2x) * sqrt(T(2))
