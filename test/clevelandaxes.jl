f = clevelandaxes!(Figure(), ["S$(lpad(i, 2))" for i in 1:16], (4, 4))
n = 12
# Hack to determine whether we're using Makie 0.19 by looking for a characteristic
# breaking change. This can be safely removed once the package requires Makie v0.19
# at a minimum.
local text!
try
    Makie.text!(Axis(Figure()[1, 1]), "test"; textsize=69)
    text! = (args...; fontsize, kwargs...) -> Makie.text!(args...;
                                                          textsize=fontsize,
                                                          kwargs...)
catch err
    if err isa ArgumentError && occursin("Makie v0.19", sprint(showerror, err))
        text! = Makie.text!
    else
        rethrow(err)
    end
end
for i in 1:4, j in 1:4
    x = randn(MersenneTwister(i), n)
    y = randn(MersenneTwister(j), n)
    scatter!(f[i, j], x, y)
    text!(f[i, j], 1.9, -1.9;
          text="[$i, $j]", align=(:center, :center), fontsize=14)
end
@test save(joinpath(OUTDIR, "clevelandaxes.png"), f)
