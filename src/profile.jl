function _mksplines(tbl)
    return map(sort(unique(tbl.i))) do i
        tbli = filter(r -> r.i == i, tbl)
        ζ = tbli.ζ
        issorted(ζ) || ArgumentError("ζ values not sorted for β[$i]")
        β = [β[i] for β in tbli.β]
        issorted(β) || ArgumentError("β values not sorted for β[$i]")
        (;
            fwd = interpolate(BSplineBasis(4, β), β, ζ),
            rev = interpolate(BSplineBasis(4, ζ), ζ, β),
        )
    end
end

function splinepanel(splpair; npts::Integer=100 )
    (; fwd) = splpair
    extr = extrema(fwd.basis.breakpoints)
    betas = collect(range(first(extr), last(extr); length=npts))
    zetas = [fwd(β) for β in betas]
    return lines(betas, zetas)
end

function zetaplot(
    pr::MixedModelProfile;
    absv::Bool=true,
    coverage=[0.5, 0.8, 0.9, 0.95, 0.99],
)
    (; prtbl, fecnames, facnames, recnames) = pr
    spl = _mksplines(prtbl)
    splinepanel(first(spl))
end
