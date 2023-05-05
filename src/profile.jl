"""
    zetaplot!(f::FigureLike, pr::MixedModelProfile;
              absv=false,
              ptyp='β',
              coverage=[.5,.8,.9,.95,.99],
              zbd=nothing)

Add axes with plots of the profile ζ (or its absolute value) for parameters
starting with `ptyp` from `pr` to `f`.

Valid `ptyp` values are 'β', 'σ', and 'θ'.

If `absv` is `true` then intervals corresponding to coverage levels in
`coverage` are added to each panel.
"""
function zetaplot!(
    f::Makie.FigureLike,
    pr::MixedModelProfile;
    absv::Bool=false,   # plot abs(zeta) vs parameter value and add intervals
    ptyp::Char='β',
    coverage=[0.5, 0.8, 0.9, 0.95, 0.99],
    zbd::Union{Nothing,Number}=nothing,
)
    ptyp in Set(['β', 'σ', 'θ']) ||
        throw(ArgumentError("Invalid `ptyp`: $(ptyp)."))
    axs = Axis[]
    cutoffs = sqrt.(quantile(Chisq(1), coverage))
    zbd = something(zbd, 1.05 * maximum(cutoffs))
    ylabel = absv ? "|ζ|" : "ζ"
    fwd, rev = pr.fwd, pr.rev
    filt = startswith(ptyp) ∘ string          # filter function
    for (i, p) in enumerate(sort!(filter(filt, collect(keys(fwd)))))
        fw, rp = fwd[p], rev[p]
        ax = Axis(f[1, i]; xlabel=string(p), ylabel)
        isone(i) || hideydecorations!(ax; grid=false)
        push!(axs, ax)
        knts = knots(rp)
        intvl = rp(max(-zbd, first(knts))) .. rp(min(zbd, last(knts)))
        lines!(ax, intvl, (absv ? abs : identity) ∘ fw)
        if absv
            sgncp = repeat(cutoffs; inner=2) .* repeat([-1, 1]; outer=length(cutoffs))
            linesegments!(
                rp.(clamp.(sgncp, first(knts), last(knts))),
                abs.(sgncp),
            )
        else
            xv = filter(x -> x ∈ intvl, fw.x)
            scatter!(ax, xv, fw.(xv))
            est = rp(0)
            slope = (Derivative(1) * fw)(est)
            ablines!(ax, -(slope * est), slope)
        end
    end
    linkyaxes!(axs...)
    f
end

"""
    zetaplot(args...; kwargs...)

Convenience wrapper for `zetaplot!(Figure(), ...)`.

See [`zetaplot!`](@ref).
"""
zetaplot(args...; kwargs...) = zetaplot!(Figure(), args...; kwargs...)

"""
    profiledensity!(f::FigureLike, pr::MixedModelProfile;
                    ptyp::Char='σ',
                    zbd=3,
                    share_y_scale=true).

Add axes with density plots of the profile ζ for parameters
starting with `ptyp` from `pr` to `f`.

Valid `ptyp` values are 'β', 'σ', and 'θ'.

If `share_y_scale`, the each facet shares a common y-scale.
"""
function profiledensity!(
    f::Makie.FigureLike,
    pr::MixedModelProfile;
    zbd=3,
    ptyp::Char='σ',
    share_y_scale=true)

    ptyp in Set(['β', 'σ', 'θ']) ||
        throw(ArgumentError("Invalid `ptyp`: $(ptyp)."))

    fwd, rev = pr.fwd, pr.rev
    ks = sort!(collect(filter(k -> startswith(string(k), ptyp), keys(fwd))))
    axs = sizehint!(Axis[], length(ks))
    for (i, p) in enumerate(ks)
        rp, fw = rev[p], fwd[p]

        ax = Axis(f[1, i]; xlabel=string(p), ylabel="pdf")
        if share_y_scale && i > 1
            hideydecorations!(ax; grid=false, ticks=false)
        end
        push!(axs, ax)
        knts = knots(rp)
        lines!(
            ax,
            rp(max(-zbd, first(knts))) .. rp(min(zbd, last(knts))),
            x -> pdf(Normal(), fw(x)) * (Derivative(1) * fw)(x)
        )
    end
    if share_y_scale
        linkyaxes!(axs...)
    end
    f
end

"""
    profiledensity(args...; kwargs...)

Convenience wrapper for `profiledensity!(Figure(), ...)`.

See [`profiledensity!`](@ref).
"""
profiledensity(args...; kwargs...) = profiledensity!(Figure(), args...; kwargs...)

#= outdated code
function zetatraces!(ax::Axis, pr::MixedModelProfile, i, j)
    (; prtbl, δ, fwd) = pr
    βmat = reshape(prtbl.β, :, length(fwd))
    i > j || throw(ArgumentError("i = $i is not greater than j = $j"))
    lines!(ax, fwd[i].(getindex.(view(βmat, :, j), i)), δ)
    lines!(ax, δ, fwd[j].(getindex.(view(βmat, :, i), j)))
    return ax
end

function zetatraceplot(pr::MixedModelProfile; figure=(; resolution=(800,800)))
    (; fecnames) = pr
    f = Figure(; figure)
    p = length(fecnames)
    for j in 1:(p - 1)
        for i in (j + 1):p
            zetatraces!(Axis(f[i - 1, j]), pr, i, j)
        end
    end
    return f
end

"""
    _signedavgdelta(x, y)

Return the signed average and half-difference of acos(x) / π and acos(y) / π

Because of the range of
"""
function _signedavgdelta(x, y)
    acpix, acpiy = acos.((x, y)) ./ π
    d = (acpix - acpiy) / 2
    return sign(d) * (acpix + acpiy) / 2, abs(d)
end

"""
    zetaijspl(pr, i, j)

Return an interpolation spline for the βj trace on βi in the zeta scale.
"""
function zetaijspline!(pr::MixedModelProfile, i, j)
    (; prtbl, δ, fwd, rev) = pr
    iinds = axes(δ, 1) .+ (i - 1) * length(δ)  # row indices for i'th trace in prtbl
    zetavalsj = fwd[j].(getindex.(view(prtbl.β, iinds), j))
    return interpolate!(rev[i], zetavalsj)
end

=#
