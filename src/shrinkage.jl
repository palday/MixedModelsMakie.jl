"""
    shrinkage(m)

Return a column table of the levels of a grouping factor, the conditional modes
of the random effects and the within-group estimates for these coefficients
"""
function shrinkage(m::LinearMixedModel{T}) where {T}
    isone(length(m.reterms)) || throw(ArgumentError("m must have a single random-effects term"))
    re = only(m.reterms)
    cnms = re.cnames
    fenms = m.feterm.cnames
    cnms ⊆ fenms || throw(ArgumentError("re.cnames is not a subset of m.feterm.cnames"))
    xcols = [findfirst(==(nm), m.feterm.cnames) for nm in cnms]
    Xmat = view(m.Xymat, :, xcols)
    yvec = view(m.Xymat, :, size(m.Xymat, 2))
    Ttyp = NamedTuple{(Symbol.(cnms)..., ), NTuple{length(cnms), T}}
    effects = only(ranef(m))
    mixed = [Ttyp(v .+ view(coef(m), xcols)) for v in eachcol(only(ranef(m)))]
    within = Vector{Union{Missing,Ttyp}}(missing, length(mixed))
    levs = re.levels
    refs = re.refs
    for i in eachindex(levs)
        rows = findall(==(i), refs)
        try
            within[i] = Ttyp(view(Xmat, rows, :) \ view(yvec, rows))
        catch
        end
    end
    (levels = levs, within = (any(ismissing.(within)) ? within : disallowmissing(within)), mixed = mixed)
end