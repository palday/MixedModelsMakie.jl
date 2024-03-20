pr1 = profile(m1)
for ptyp in ['σ', 'θ', 'β'], toggle in [true, false]
    local f
    f = zetaplot(pr1; ptyp, absv=toggle)
    save(joinpath(OUTDIR, "zetaplot_$(ptyp)_$(toggle).png"), f)

    f = profiledensity(pr1; ptyp, share_y_scale=toggle)
    save(joinpath(OUTDIR, "profiledensity_$(ptyp)_$(toggle).png"), f)
end
