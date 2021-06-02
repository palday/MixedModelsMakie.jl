"""
    clevelandaxes!(f::Figure, labs, layout)

Create a set of axes within `f` with rows, columns determined by `layout` sufficient to hold `labs`
"""
function clevelandaxes!(f::Figure, labs, layout)
    nrow, ncol = layout
    npanel = length(labs)
    axs = sizehint!(Axis[], npanel)
    for ind in eachindex(labs)
        i, j = fldmod1(ind, ncol)
        ii = nrow + 1 - i   # axes grid rows numbered from top, we want to fill bottom row first
        ax = Axis(f[ii, j])
        push!(axs, ax)
        if isone(i)         # on the bottom row
            hidexdecorations!(ax, grid=false, ticks=false, ticklabels=iseven(j))
        elseif isone(ii)    # on the top row
            hidexdecorations!(ax, grid=false, ticks=false, ticklabels=isodd(j))
            ax.xaxisposition = :top
        else
            hidexdecorations!(ax, grid=false)
        end
        if isone(j)         # on the left side
            hideydecorations!(ax, grid=false, ticks=false, ticklabels=iseven(i))
        elseif j == ncol    # on the right side
            hideydecorations!(ax, grid=false, ticks=false, ticklabels=isodd(i))
            ax.yaxisposition = :right
        else
            hideydecorations!(ax, grid=false)
        end
    end
    linkaxes!(axs...)
    colgap!(f.layout, 0)
    rowgap!(f.layout, 0)
    f
end

"""
    simplelinreg(x, y)

Return a Tuple of the coefficients, `(a, b)`,  from a simple linear regression, `y = a + bx + ϵ`
"""
function simplelinreg(x, y)
    x, y = float(x), float(y)
    A = cholesky!(Symmetric([length(x) sum(x) sum(y); 0.0 sum(abs2, x) dot(x, y); 0.0 0.0 sum(abs2, y)])).factors
    (ldiv!(UpperTriangular(view(A, 1:2, 1:2)), view(A, 1:2, 3))..., )
end
