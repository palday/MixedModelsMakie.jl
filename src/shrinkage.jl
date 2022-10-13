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
        return linkxaxes!(col...)
    end

    return f
end

getellipsepoints(radius, lambda) = getellipsepoints(0, 0, radius, lambda)

function getellipsepoints(cx, cy, radius, lambda)
    t = range(0, 2π; length=100)
    ellipse_x_r = cos.(t)
    ellipse_y_r = sin.(t)
    r_ellipse = radius .* hcat(ellipse_x_r, ellipse_y_r) * lambda'
    x = @. cx + r_ellipse[:, 2]
    y = @. cy + r_ellipse[:, 1]
    return x, y
end

"""
    shrinkageplot!(f::Figure, m::MixedModel, gf::Symbol=first(fnames(m)), θref;
                   ellipse=false, ellipse_scale=1, n_ellipse=5)

Return a scatter-plot matrix of the conditional means, b, of the random effects for grouping factor `gf`.

Two sets of conditional means are plotted: those at the estimated parameter values and those at `θref`.
The default `θref` results in `Λ` being a very large multiple of the identity.  The corresponding
conditional means can be regarded as unpenalized.

Correlation ellipses can be added with `ellipse=true`, with the number of ellipses controlled by
`n_ellipse`. The ellipses are equally spaced between the outer ellipse and the origin (center).
The scaling of the ellipses can be adjusted with the multiplicative `ellipse_scale`. If you are
unable to see the ellipses, try increasing `ellipse_scale`.

!!! note
    For degenerate (singular) models, the correlation ellipse will also be degenerate, i.e.,
    collapse to a point or line.
"""
function shrinkageplot!(f::Figure,
                        m::MixedModel{T},
                        gf::Symbol=first(fnames(m)),
                        θref::AbstractVector{T}=(isa(m, LinearMixedModel) ? 1e4 : 1) .*
                                                m.optsum.initial;
                        ellipse=false, ellipse_scale=1, n_ellipse=5) where {T}
    reind = findfirst(==(gf), fnames(m))  # convert the symbol gf to an index
    if isnothing(reind)
        throw(ArgumentError("gf=$gf is not one of the grouping factor names, $(fnames(m))"))
    end
    reest = ranef(m)[reind]          # random effects conditional means at estimated θ
    reref = _ranef(m, θref)[reind]   # same at θref
    remat = m.reterms[reind]
    # display(remat.λ)
    function pfunc(ax, i, j)
        x, y = view(reref, j, :), view(reref, i, :)
        u, v = view(reest, j, :), view(reest, i, :)
        if ellipse
            cho = remat.λ[[i, j], [j, i]]
            rad_outer = ellipse_scale * max(maximum(abs, x), maximum(abs, y))
            rad_inner = 0
            for radius in LinRange(rad_inner, rad_outer, n_ellipse + 1)
                ex, ey = getellipsepoints(radius, cho)
                lines!(ax, ex, ey; color=:green, linestyle=:dash)
            end
        end
        scatter!(ax, x, y; color=(:red, 0.25))   # reference points
        arrows!(ax, x, y, u .- x, v .- y)        # first so arrow heads don't obscure pts
        return scatter!(ax, u, v; color=(:blue, 0.25))  # conditional means at estimates
    end
    splomaxes!(f, m.reterms[reind].cnames, pfunc)

    return f
end

function _ranef(m::LinearMixedModel, θref)
    vv = try
        ranef(updateL!(setθ!(m, θref)))
    catch e
        @error "Failed to compute unshrunken values with the following exception:"
        rethrow(e)
    finally
        updateL!(setθ!(m, m.optsum.final)) # restore parameter estimates and update m
    end
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
    shrinkageplot(m::MixedModel, gf::Symbol=first(fnames(m)), θref, args...; kwargs...)

Return a scatter-plot matrix of the conditional means, b, of the random effects for grouping factor `gf`.

Two sets of conditional means are plotted: those at the estimated parameter values and those at `θref`.
The default `θref` results in `Λ` being a very large multiple of the identity.  The corresponding
conditional means can be regarded as unpenalized.

`args...` and `kwargs...` are passed on to [`shrinkageplot!`](@ref)
"""
function shrinkageplot(m::MixedModel, args...; kwargs...)
    f = Figure(; resolution=(1000, 1000)) # use an aspect ratio of 1 for the whole figure

    return shrinkageplot!(f, m, args...; kwargs...)
end

"""
    splom!(f::Figure, df::DataFrame)

Create a scatter-plot matrix in `f` from the columns of `df`.

Non-numeric columns are ignored.
"""
function splom!(f::Figure, df::DataFrame; addcontours::Bool=false)
    n_cols = ncol(df)
    df = select(df, findall(col -> eltype(col) <: Number, eachcol(df));
                copycols=false)
    n_cols > ncol(df) &&
        @info "Ignoring $(n_cols - ncol(df)) non-numeric columns."
    mat = Array(df)
    function pfunc(ax, i, j)
        v = view(mat, :, [j, i])
        scatter!(ax, v; color=(:blue, 0.2))
        addcontours && contour!(ax, kde(v))
        return ax
    end
    return splomaxes!(f, names(df), pfunc)
end
