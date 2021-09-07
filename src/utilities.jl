"""
    confint_table(x::StatsBase.StatisticalModel, level=0.95)
    confint_table(x::MixedModelBootstrap, level=0.95)

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

"""
function confint_table(x::StatsBase.StatisticalModel, level=0.95)
    # taking from the lower tail
    semultiple = zquantile((1 - level) / 2)
    se = stderror(x)

    return DataFrame(;
        coefname=coefnames(x),
        estimate=coef(x),
        lower=coef(x) + semultiple * se,
        upper=coef(x) - semultiple * se,
    )
end

function confint_table(x::MixedModelBootstrap, level=0.95)
    df = transform!(
        select!(DataFrame(x.β), Not(:iter)), :coefname => ByRow(string) => :coefname
    )
    return combine(
        groupby(df, :coefname),
        :β => mean => :estimate,
        :β => NamedTuple{(:lower, :upper)} ∘ shortestcovint => [:lower, :upper],
    )
end


"""
    ppoints(n::Integer)

Return a sequence of `n` equally-spaced points in the interval (0, 1) - so-called "probability points"
"""
ppoints(n::Integer) = inv(2n):inv(n):1

"""
    zquantile(x::AbstractFloat)

Evaluate `quantile(Normal(), x)` using only the `SpecialFunctions` package (i.e. not requiring `Distributions`).
"""
zquantile(x::T) where {T<:AbstractFloat} = -erfcinv(2x) * sqrt(T(2))
