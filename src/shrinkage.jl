"""
    CoefByGroup

Coefficients for groups of observations defined by levels of a grouping factor

Fields include:
- `cnames`: `Vector{String}` of the column (or coefficient) names
- `levels`: `Vector` of the levels of the grouping factor.  Usually `String`s but not necessarily.
- `fixed`: `AbstractVector{T}` of the global OLS estimates for the coefficients in `cnames` only.
- `condmodes`: `AbstractMatrix{T}` of size `length(levels) x length(cnames)` of the conditional means/modes for the random effects
- `grpest`: similar to `condmodes` but allowing for `Missing` values, giving the within group OLS estimates of the coefficients.
"""
struct CoefByGroup{T<:AbstractFloat}
    cnames::Vector{String}
    levels::Vector
    fixed::AbstractVector{T}
    condmodes::AbstractMatrix{T}
    grpest::AbstractMatrix{Union{Missing,T}}
end

"""
    shrinkage(m)

Return a `NamedTuple{fnames(m), NTuple(k, CoefByGroup)}` of the coefficients by group for the grouping factors
"""
function shrinkage(m::LinearMixedModel{T}) where {T}
    fenms = m.feterm.cnames
    yvec = view(m.Xymat, :, size(m.Xymat, 2))
    ranefs = ranef(m)
    cvec = sizehint!(CoefByGroup[], length(ranefs))
    for (j, re) in enumerate(m.reterms)
        cnms = re.cnames
        cnms ⊆ fenms || throw(ArgumentError("re[$j].cnames is not a subset of m.feterm.cnames"))
        Xmat = view(m.Xymat, :, [findfirst(==(nm), fenms) for nm in cnms])
        levs = re.levels
        refs = re.refs
        grpest = Matrix{Union{Missing,T}}(missing, (length(levs), length(cnms)))
        for i in eachindex(levs)
            rows = findall(==(i), refs)
            try
                grpest[i, :] = view(Xmat, rows, :) \ view(yvec, rows)
            catch
            end
        end
        push!(cvec, CoefByGroup(cnms, levs, Xmat \ yvec, ranefs[j]', grpest))
    end
    NamedTuple{fnames(m)}((cvec...,))
end

shrinkageplot(m::LinearMixedModel, gf::Symbol=first(fnames(m))) = shrinkageplot(shrinkage(m)[gf])

shrinkageplot(cg::CoefByGroup) = shrinkageplot(Figure(resolution=(1000,1000)), cg)

function shrinkageplot(f::Figure, cg::CoefByGroup)
    shrinkage2d(Axis(f[1,1]), cg)
    f
end

function shrinkage2d(a::Axis, cg::CoefByGroup{T}, inds=(1, 2)) where T
    i, j = inds
    fxd = cg.fixed
    scatter!(a, view(fxd, i:i), view(fxd, j:j), color=:green, label="Pop")
    u = view(cg.condmodes, :, i) .+ fxd[i]
    v = view(cg.condmodes, :, j) .+ fxd[j]
    scatter!(a, u, v, color=:blue, label="Mixed")
    x = view(cg.grpest, :, i)
    y = view(cg.grpest, :, j)
    missgrp = ismissing.(x) .| ismissing.(y)
    if any(missgrp)
        nonmiss = .!(missgrp)
        x = convert(Vector{T}, view(x, nonmiss))
        y = convert(Vector{T}, view(y, nonmiss))
        u = view(u, nonmiss)
        v = view(v, nonmiss)
    else
        x = convert(Vector{T}, x)
        y = convert(Vector{T}, y)
    end
    scatter!(a, x, y, color=:red, label="Within")
    arrows!(a, x, y, u - x, v - y)
    a.xlabel = cg.cnames[i]
    a.ylabel = cg.cnames[j]
    a
end
