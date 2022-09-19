function zetaplot(
    pr::MixedModelProfile;
    absv::Bool=true,
    coverage=[0.5, 0.8, 0.9, 0.95, 0.99],
)
    (; prtbl, δ, fecnames) = pr
    draw(
        data(
            (;
            δ = repeat(collect(δ); outer=length(fecnames)),
            ζ = prtbl.ζ,
            cnames = repeat(fecnames, inner=length(δ)),
            ),
        ) *
        (
            if absv
                mapping(:δ, :ζ => abs => "|ζ|"; color=(:cnames => "Coefficient"))
            else
                mapping(:δ, :ζ; color=(:cnames => "Coefficient"))
            end
         ) *
        visual(Lines);
        figure=(; resolution=(1000, 600)),
    )
end

function deltadensity(pr::MixedModelProfile)
    (; δ, fecnames, splines) = pr
    xpat = first(δ):0.125:last(δ)
    xv = repeat(xpat, outer=length(fecnames))
    dv = @. exp(-abs2(xv) / 2) / inv(sqrt(2π)) 
    dens = dv ./ foldl(vcat, map(s -> (Derivative(1) * s).(xpat), splines))
    cnames = repeat(fecnames; inner=length(xpat))
    draw( 
        data((; xv, dens, cnames)) *
        mapping(
            :xv => "δ",
            :dens => "Effective probability density";
            color=:cnames => "Coefficient",
        ) * 
        visual(Lines);
        figure=(; resolution=(1000, 600)),
    )
end

