function zetaplot(
    pr::MixedModelProfile;
    absv::Bool=true,
    coverage=[0.5, 0.8, 0.9, 0.95, 0.99],
    resolution=(1200,600),
)
    cutoffs = sqrt.(quantile(Chisq(1), coverage))
    zbd = 1.05 * maximum(cutoffs)
    f = Figure(; resolution)
    fwd, rev = pr.fwd, pr.rev
    ylabel = absv ? "|ζ|" : "ζ"
    for (i, p) in enumerate(keys(fwd))
        rp, fw = rev[p], fwd[p]
        ax = Axis(f[1, i]; xlabel=string(p), ylabel)
        knts = knots(rp)
        intvl = rp(max(-zbd, first(knts)))..rp(min(zbd, last(knts)))
        lines!(ax, intvl, (absv ? abs : identity) ∘ fw)
        if absv
            sgncp = repeat(cutoffs; inner=2) .* repeat([-1, 1]; outer=length(cutoffs))
            linesegments!(
                rp.(clamp.(sgncp, first(knts), last(knts))),
                abs.(sgncp),
            )
        end
    end
    f
end

function profiledensity(
    pr::MixedModelProfile;
    resolution=(1200, 500),
    zbd=3,
)
    fwd, rev = pr.fwd, pr.rev
    f = Figure(; resolution)
    for (i, p) in enumerate(keys(fwd))
        rp, fw = rev[p], fwd[p]
        ax = Axis(f[1, i]; xlabel=string(p), ylabel="pdf")
        knts = knots(rp)
        intvl = rp(max(-zbd, first(knts)))..rp(min(zbd, last(knts)))
        lines!(
            ax,
            intvl,
            x -> pdf(Normal(), fw(x)) * (Derivative(1) * fw)(x) 
        )
    end
    f
end

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
