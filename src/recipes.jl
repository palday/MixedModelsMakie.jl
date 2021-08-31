# this file is for short recipes

# XXX it would be great to have a 1-1 aspect ratio here,
# but this seems like something that should be done upstream
function Makie.convert_arguments(P::Type{<:Makie.QQNorm}, x::MixedModel, args...; kwargs...)
    return convert_arguments(P, residuals(x), args...; kwargs...)
end
function Makie.convert_arguments(
    P::Type{<:Makie.QQPlot}, d::Distributions.Distribution, x::MixedModel, args...; kwargs...
)
    return convert_arguments(P, d, residuals(x), args...; kwargs...)
end

# Makie.convert_arguments(P::PointBased, x::MixedModel) = convert_arguments(P, response(x), fitted(x))
# Makie.convert_arguments(P::Type{<:Makie.Scatter}, x::MixedModel) = convert_arguments(P, fitted(x), residuals(x))
# Makie.plottype(::MixedModel) = Scatter
