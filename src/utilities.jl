function _coefnames(x; show_intercept=true)
    cn = fixefnames(x)
    return show_intercept ? cn : filter!(!=("(Intercept)"), cn)
end

function _coefnames(x::MixedModelBootstrap, ptype; show_intercept=true)
    nt = getproperty(first(x.fits), ptype)
    cn = [string(k) for (k, v) in pairs(nt) if !isequal(v, -0.0)]
    return show_intercept ? cn : filter!(!=("(Intercept)"), cn)
end

"""
    confint_table(x::MixedModel, level=0.95; show_intercept=true)
    confint_table(x::MixedModelBootstrap, level=0.95; show_intercept=true)

Return a DataFrame of coefficient names, point estimates and confidence intervals.

`level` specifies the confidence level, e.g. `0.95` for 95% confidence intervals.

For `MixedModels`, the intervals are computed using the standard errors and the asymptotic
Wald approximation (e.g. est±1.96*se for 95% intervals). For `MixedModelBootstrap`, the intervals
are computed using `shortestcovint`, and the point estimate re-computed as the mean of the bootstrap values.

The returned table has the following columns:
- `coefname`: the names of the coefficients
- `estimate`: the point estimates
- `lower`: the lower edge of the confidence interval
- `upper`: the upper edge of the confidence interval

!!! note
    This function is internal and may be removed in a future release
    without being considered a breaking change.
"""
function confint_table(x::MixedModel, level=0.95; show_intercept=true)
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

function confint_table(x::MixedModelBootstrap, level=0.95; show_intercept=true)
    df = transform!(select!(DataFrame(x.β), Not(:iter)),
                    :coefname => ByRow(string) => :coefname)
    ci(x) = shortestcovint(x, level)
    df = combine(groupby(df, :coefname),
                 :β => mean => :estimate,
                 :β => NamedTuple{(:lower, :upper)} ∘ ci => [:lower, :upper])
    return filter!(:coefname => in(_coefnames(x; show_intercept)), df)
end

_npreds(x; show_intercept=true) = length(_coefnames(x; show_intercept))

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
