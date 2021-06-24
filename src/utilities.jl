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
