"""
    clevelandaxes!(f::Union{Makie.FigureLike,Makie.GridLayout}, labs, layout)

Create a set of axes within `f` with rows, columns determined by `layout` sufficient to hold `labs`
"""
function clevelandaxes!(f::Union{Makie.FigureLike,Makie.GridLayout}, labs, layout)
    nrow, ncol = layout
    npanel = length(labs)
    axs = sizehint!(Axis[], npanel)
    for ind in eachindex(labs)
        i, j = fldmod1(ind, ncol)
        ii = nrow + 1 - i   # axes grid rows numbered from top, we want to fill bottom row first
        ax = Axis(f[ii, j])
        push!(axs, ax)
        if isone(i)         # on the bottom row
            hidexdecorations!(ax; grid=false, ticks=false, ticklabels=iseven(j))
        elseif isone(ii)    # on the top row
            hidexdecorations!(ax; grid=false, ticks=false, ticklabels=isodd(j))
            ax.xaxisposition = :top
        else
            hidexdecorations!(ax; grid=false)
        end
        if isone(j)         # on the left side
            hideydecorations!(ax; grid=false, ticks=false, ticklabels=iseven(i))
        elseif j == ncol    # on the right side
            hideydecorations!(ax; grid=false, ticks=false, ticklabels=isodd(i))
            ax.yaxisposition = :right
        else
            hideydecorations!(ax; grid=false)
        end
    end
    linkaxes!(axs...)
    colgap!(f.layout, 0)
    rowgap!(f.layout, 0)
    return f
end

"""
    simplelinreg(x, y)

Return a Tuple of the coefficients, `(a, b)`,  from a simple linear regression, `y = a + bx + ϵ`
"""
function simplelinreg(x, y)
    x, y = float(x), float(y)
    A = cholesky!(Symmetric([length(x) sum(x) sum(y); 0.0 sum(abs2, x) dot(x, y);
                             0.0 0.0 sum(abs2, y)])).factors
    return (ldiv!(UpperTriangular(view(A, 1:2, 1:2)), view(A, 1:2, 3))...,)
end

"""
    splom!(f::Union{Makie.FigureLike,Makie.GridLayout}, df::DataFrame)

Create a scatter-plot matrix in `f` from the columns of `df`.

Non-numeric columns are ignored.
"""
function splom!(f::Union{Makie.FigureLike,Makie.GridLayout}, df::DataFrame;
                addcontours::Bool=false)
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

"""
    splomaxes!(f::Union{Makie.FigureLike,Makie.GridLayout}, labels::AbstractVector{<:AbstractString},
               panel!::Function, args...;
               extraticks::Bool=false, kwargs...)

Populate f with a set of `(k*(k-1))/2` axes in a lower triangle for all pairs of `labels`,
where `k` is the length of `labels`.  The `panel!` function should have the signature
`panel!(ax::Axis, i::Integer, j::Integer, args...; kwargs...)` and should draw the
[i,j] panel in `ax`.
"""
function splomaxes!(f::Union{Makie.FigureLike,Makie.GridLayout},
                    labels::AbstractVector{<:AbstractString},
                    panel!::Function, args...; extraticks::Bool=false, kwargs...)
    k = length(labels)
    cols = Dict()
    for i in 2:k                          # strict lower triangle of panels
        row = Axis[]
        for j in 1:(i - 1)
            ax = Axis(f[i - 1, j])
            panel!(ax, i, j, args...; kwargs...)
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
