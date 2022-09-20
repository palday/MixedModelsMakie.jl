function zetaplot(
    pr::MixedModelProfile;
    absv::Bool=true,
    coverage=[0.5, 0.8, 0.9, 0.95, 0.99],
    figure=(; resolution=(1000,600))
)
    cutoffs = sqrt.(quantile(Chisq(1), coverage))
    yaxs = absv ? :ζ => abs => "|ζ|" : :ζ
    (; prtbl, δ, fecnames) = pr
    draw(
        data(
            (;
            δ = repeat(collect(δ); outer=length(fecnames)),
            ζ = prtbl.ζ,
            cnames = repeat(fecnames, inner=length(δ)),
            ),
        ) *
        mapping(:δ, yaxs, color=:cnames => "Coefficient") *
        visual(Lines);
        figure,
    )
end

function deltadensity(
    pr::MixedModelProfile;
    npts::Integer=129,
    figure=(; resolution=(1000, 600)),
)
    (; prtbl, δ, fecnames) = pr
    ord, nat = BSplineOrder(4), Natural()
    splines = [interpolate(ζ, δ, ord, nat) for ζ in eachcol(reshape(prtbl.ζ, length(δ), :))]
    xpat = range(first(δ), last(δ); length=npts)
    xv = repeat(xpat; outer=length(fecnames))
    cnames = repeat(fecnames; inner=length(xpat))
    dv = @. exp(-abs2(xv) / 2) / inv(sqrt(2π)) 
    dens = dv .* foldl(vcat, map(s -> (Derivative(1) * s).(xpat), splines))
    draw( 
        data((; xv, dens, cnames)) *
        mapping(
            :xv => "δ",
            :dens => "Effective probability density";
            color=:cnames => "Coefficient",
        ) * 
        visual(Lines);
        figure,
    )
end

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

