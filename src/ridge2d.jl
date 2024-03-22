function _ridge2d_panel!(ax::Axis, i::Int, j::Int, cnames::Vector{Symbol}, tbl)
    x = Tables.getcolumn(tbl, cnames[j])
    y = Tables.getcolumn(tbl, cnames[i])
    dens = kde((x, y))
    scatter!(ax, x, y; color=:black, alpha=0.2)
    plt = contour!(ax, collect(dens.x), collect(dens.y), dens.density;
                   color=:green, linewidth=3,
                   labelsize=30, labels=false)   # reference points

    return plt
end

"""
    ridge2d!(f::Union{Makie.FigureLike,Makie.GridLayout}, bs::MixedModelBootstrap; 
             ptype=:β)

Plot pairwise bivariate scatter plots with overlain densities for a bootstrap sample.


`ptype` specifies the set of parameters to examine, e.g. `:β`, `:σ`, `:ρ`.
"""
function ridge2d!(f::Indexable, bs::MixedModelBootstrap;
                  ptype=:β)
    tbl = bs.tbl
    cnames = [string(x) for x in propertynames(tbl)[2:end]]
    filter!(startswith(string(ptype)), cnames)
    isempty(cnames) &&
        throw(ArgumentError("No parameters $ptype found."))
    length(cnames) == 1 &&
        throw(ArgumentError("Only 1 $ptype-parameter found: 2D plots require at least 2."))
    splomaxes!(f, cnames, _ridge2d_panel!, Symbol.(cnames), tbl)
    return f
end

"""$(@doc ridge2d!)"""
function ridge2d(bs::MixedModelBootstrap, args...; kwargs...)
    f = Figure(; size=(1000, 1000)) # use an aspect ratio of 1 for the whole figure

    return ridge2d!(f, bs, args...; kwargs...)
end
