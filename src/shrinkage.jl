"""
    CoefByGroup

Coefficients for groups of observations defined by levels of a grouping factor

Fields include:
- `globalest`: `NamedTuple` - names are coefficient names, values are the global estimates
- `withinmxdtbl`: `Tables.columntable` (i.e. `NamedTuple` of `Vectors`).  Names are `level`, `mixed` and `within`

!!! note
    This functionality may be moved upstream into MixedModels.jl in the near future.
"""
struct CoefByGroup
    globalest::NamedTuple    # global coefficient estimates
    withinmxdtbl::NamedTuple # Tables.columntable of levels, within-group estimates, and mixed-effects "estimates"
end

"""
    shrinkage(m)

Return a `NamedTuple{fnames(m), NTuple(k, CoefByGroup)}` of the coefficients by group for the grouping factors
"""
function shrinkage(m::LinearMixedModel{T}) where {T}
    # should this be fixefnames or coefnames?
    # need to think about when pivoting is occuring
    fenms = fixefnames(m)
    # returns a view from Xymat in MM4, or just the response vec in MM3.x
    yvec = response(m)
    ranefs = ranef(m)
    cvec = []
    for (j, re) in enumerate(m.reterms)
        cnms, levs, refs = re.cnames, re.levels, re.refs
        cnms ⊆ fenms || throw(ArgumentError("No corresponding fixed effect for random effects $(setdiff(cnms, fenms)): there is no estimated grand mean to measure shrinkage towards"))
        # XXX Should m.X return a view into m.Xymat?  It could but needs MM4 to do so
        # Is it a big deal to take a view of a view?
        Xmat = view(m.X, :, [findfirst(==(nm), fenms) for nm in cnms])
        globalvec = Xmat \ yvec
        globalest = NamedTuple{(Symbol.(cnms)...,)}((globalvec..., ))
        memat = ranefs[j] .+ globalvec
        meest = [(view(memat, :, i)...,) for i in axes(memat, 2)]
        grpest = NTuple{length(cnms), T}[]
        for i in eachindex(levs)
            rows = findall(==(i), refs)
            push!(grpest, ((view(Xmat, rows, :) \ view(yvec, rows))...,))
        end
        push!(cvec, CoefByGroup(globalest, (level=levs, mixed=meest, within=grpest)))
    end
    NamedTuple{fnames(m)}((cvec...,))
end

shrinkageplot(m::LinearMixedModel, gf::Symbol=first(fnames(m))) = shrinkageplot(shrinkage(m)[gf])

shrinkageplot(cg::CoefByGroup) = shrinkageplot!(Figure(resolution=(1000,1000)), cg)

function shrinkageplot!(f::Figure, cg::CoefByGroup, inds=(1,2))
    length(inds) == 2 && shrinkage2d!(Axis(f[1,1]), cg, inds)
    f
end

function shrinkage2d!(a::Axis, cg::CoefByGroup, inds=(1, 2))
    fxd = cg.globalest
    T = typeof(first(fxd))
    i, j = inds
    scatter!(a, [fxd[i]], [fxd[j]], color=:green, label="Pop")
    wmx = cg.withinmxdtbl
    u = getindex.(wmx.mixed, i)
    v = getindex.(wmx.mixed, j)
    scatter!(a, u, v, color=:blue, label="Mixed")
    x = getindex.(wmx.within, i)
    y = getindex.(wmx.within, j)
    scatter!(a, x, y, color=:red, label="Within")
    arrows!(a, x, y, u - x, v - y)
    nms = keys(cg.globalest)
    a.xlabel = string(nms[i])
    a.ylabel = string(nms[j])
    a
end
