m6 = let
    dat = MixedModels.dataset(:mrk17_exp1)
    form = @formula(1000 / rt ~ 1 + F*P*Q*lQ*lT + (1+F+P+Q+lQ+lT|subj) + (1+P+Q+lQ+lT|item))
    contr = Dict(
        :F => EffectsCoding(),
        :P => EffectsCoding(),
        :Q => EffectsCoding(),
        :lQ => EffectsCoding(),
        :lT => EffectsCoding(),
        :subj => Grouping(),
        :item => Grouping(),
    )
    fit(MixedModel, form, dat; contrasts = contr)
end
