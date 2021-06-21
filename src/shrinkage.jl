"""
    shrinkageplot(m::LinearMixedModel, gf::Symbol=first(fnames(m)), θref)

Return a scatter-plot matrix of the conditional means, b, of the random effects for grouping factor `gf`.

Two sets of conditional means are plotted: those at the estimated parameter values and those at `θref`.
The default `θref` results in `Λ` being a very large multiple of the identity.  The corresponding
conditional means can be regarded as unpenalized.
"""
function shrinkageplot(
    m::LinearMixedModel{T},
    gf::Symbol=first(fnames(m)),
    θref::AbstractVector{T}=10000 .* m.optsum.initial,
) where {T}
    reind = findfirst(==(gf), fnames(m))  # convert the symbol gf to an index
    if isnothing(reind)
        throw(ArgumentError("gf=$gf is not one of the grouping factor names, $(fnames(m))"))
    end
    reest = ranef(m)[reind]               # random effects conditional means at estimated θ
    reref = ranef(updateL!(setθ!(m, θref)))[reind]  # same at θref
    updateL!(setθ!(m, m.optsum.final))    # restore parameter estimates and update m
    cnms = m.reterms[reind].cnames
    f = Figure(; resolution=(1000, 1000)) # use an aspect ratio of 1 for the whole figure
    k = size(reest, 1)  # dimension of the random-effects vector per level of gf
    for i in 2:k                          # strict lower triangle of panels
        for j in 1:(i - 1)
            ax = Axis(f[i - 1, j])
            x, y = view(reref, j, :), view(reref, i, :)
            scatter!(ax, x, y; color=(:red, 0.25))   # reference points
            u, v = view(reest, j, :), view(reest, i, :)
            arrows!(ax, x, y, u .- x, v .- y)        # first so arrow heads don't obscure pts
            scatter!(ax, u, v; color=(:blue, 0.25))  # conditional means at estimates
            if i == k              # add x labels on bottom row
                ax.xlabel = string(cnms[j])
            else
                hidexdecorations!(ax; grid=false)
            end
            if isone(j)            # add y labels on left column
                ax.ylabel = string(cnms[i])
            else 
                hideydecorations!(ax; grid=false)
            end
        end
    end
    return f
end
