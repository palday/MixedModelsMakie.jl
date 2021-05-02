function Clevelandaxes!(f::Figure, npanel, layout)
    nrow, ncol = layout
    axs = sizehint!(Axis[], npanel)
    for ind in 1:npanel
        i, j = fldmod1(ind, ncol)
        ii = nrow + 1 - i   # axes grid rows numbered from top, we want to fill bottom row first
        ax = Axis(f[ii, j])
        push!(axs, ax)
        if isone(i)  # on the bottom row
            hidexdecorations!(ax, grid=false, ticks=false, ticklabels=iseven(j))
        elseif isone(ii)
            hidexdecorations!(ax, grid=false, ticks=false, ticklabels=isodd(j))
            ax.xaxisposition = :top
        else
            hidexdecorations!(ax, grid=false)
        end
        if isone(j)
            hideydecorations!(ax, grid=false, ticks=false, ticklabels=iseven(i))
        elseif j == ncol
            hideydecorations!(ax, grid=false, ticks=false, ticklabels=isodd(i))
            ax.yaxisposition = :right
        else
            hideydecorations!(ax, grid=false)
        end
    end
    linkaxes!(axs...)
    f
end

"""
    simplelinreg(x, y)

Return the coefficients, `[a, b]`,  from a simple linear regression, `y = a + bx + ϵ`
"""
function simplelinreg(x, y)
    A = cholesky!(Symmetric([length(x) sum(x) sum(y); 0.0 sum(abs2, x) dot(x, y); 0.0 0.0 sum(abs2, y)])).factors
    ldiv!(UpperTriangular(view(A, 1:2, 1:2)), view(A, 1:2, 3))
end