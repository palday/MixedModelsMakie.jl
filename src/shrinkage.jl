"""
    splomaxes!(f::Figure, labels::AbstractVector{<:AbstractString})

Populate f with a set of `(k*(k-1))/2` axes in a lower triangle for all pairs of `labels`,
where `k` is the length of `labels`.  The `panel!` function should have the signature
`panel!(ax::Axis, i::Integer, j::Integer)` and should draw the [i,j] panel in `ax`.
"""
function splomaxes!(f::Figure, labels, panel!::Function; extraticks::Bool=false)
    k = length(labels)
    cols = Dict()
    for i in 2:k                          # strict lower triangle of panels
        row = Axis[]
        for j in 1:(i - 1)
            ax = Axis(f[i - 1, j])
            panel!(ax, i, j)
            push!(row, ax)
            col = get!(cols, j, Axis[])
            push!(col, ax)
            if i == k              # add x labels on bottom row
                ax.xlabel = string(labels[j])
            elseif extraticks && i == 2
                ax.xaxisposition = :top
                hidexdecorations!(ax; grid=false, ticks=false)
            else
                hidexdecorations!(ax; grid=false)
            end
            if isone(j)            # add y labels on left column
                ax.ylabel = string(labels[i])
            elseif extraticks && j == i - 1
                ax.yaxisposition = :right
                hideydecorations!(ax; grid=false, ticks=false)
            else
                hideydecorations!(ax; grid=false)
            end
        end
        linkyaxes!(row...)
    end

    foreach(values(cols)) do col
        linkxaxes!(col...)
    end

    return f
end

"""
    shrinkageplot!(f::Figure, m::MixedModel, gf::Symbol=first(fnames(m)), θref)

Return a scatter-plot matrix of the conditional means, b, of the random effects for grouping factor `gf`.

Two sets of conditional means are plotted: those at the estimated parameter values and those at `θref`.
The default `θref` results in `Λ` being a very large multiple of the identity.  The corresponding
conditional means can be regarded as unpenalized.
"""
function shrinkageplot!(
    f::Figure,
    m::MixedModel{T},
    gf::Symbol=first(fnames(m)),
    θref::AbstractVector{T}= (isa(m, LinearMixedModel) ? 1e4 : 1) .* m.optsum.initial,
) where {T}
    reind = findfirst(==(gf), fnames(m))  # convert the symbol gf to an index
    if isnothing(reind)
        throw(ArgumentError("gf=$gf is not one of the grouping factor names, $(fnames(m))"))
    end
    reest = ranef(m)[reind]          # random effects conditional means at estimated θ
    reref = _ranef(m, θref)[reind]   # same at θref
    function pfunc(ax, i, j)
        x, y = view(reref, j, :), view(reref, i, :)
        scatter!(ax, x, y; color=(:red, 0.25))   # reference points
        u, v = view(reest, j, :), view(reest, i, :)
        arrows!(ax, x, y, u .- x, v .- y)        # first so arrow heads don't obscure pts
        return scatter!(ax, u, v; color=(:blue, 0.25))  # conditional means at estimates
    end
    splomaxes!(f, m.reterms[reind].cnames, pfunc)

    return f
end

function _ranef(m::LinearMixedModel, θref)
    vv = ranef(updateL!(setθ!(m, θref)))
    updateL!(setθ!(m, m.optsum.final)) # restore parameter estimates and update m
    return vv
end

function _ranef(m::GeneralizedLinearMixedModel, θref)
    fast = length(m.θ) == length(m.optsum.final)
    setpar! = fast ? MixedModels.setθ! : MixedModels.setβθ!
    vv = try
        ranef(pirls!(setpar!(m, θref), fast, false)) # not verbose
    catch e
        @error "Failed to compute unshrunken values with the following exception:"
        rethrow(e)
    finally
        pirls!(setpar!(m, m.optsum.final), fast, false) # restore parameter estimates and update m
    end

    return vv
end

"""
    shrinkageplot(m::MixedModel, gf::Symbol=first(fnames(m)), θref)

Return a scatter-plot matrix of the conditional means, b, of the random effects for grouping factor `gf`.

Two sets of conditional means are plotted: those at the estimated parameter values and those at `θref`.
The default `θref` results in `Λ` being a very large multiple of the identity.  The corresponding
conditional means can be regarded as unpenalized.
"""
function shrinkageplot(m::MixedModel, args...)
    f = Figure(; resolution=(1000, 1000)) # use an aspect ratio of 1 for the whole figure

    return shrinkageplot!(f, m, args...)
end

"""
    splom!(f::Figure, df::DataFrame)

Create a scatter-plot matrix in `f` from the columns of `df`.
"""
function splom!(f::Figure, df::DataFrame; addcontours::Bool=false)
    mat = Array(df)
    function pfunc(ax, i, j)
        v = view(mat, :, [j, i])
        scatter!(ax, v, color=(:blue, 0.2))
        addcontours && contour!(ax, kde(v))
        return ax
    end
    splomaxes!(f, names(df), pfunc)
end
