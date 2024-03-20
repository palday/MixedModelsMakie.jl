df = DataFrame(MixedModels.dataset(:mmec))
splof = @test_logs (:info,
                    r"Ignoring 3 non-numeric columns") splom!(Figure(), df)
@test save(joinpath(OUTDIR, "splom_mmec.png"), splof)
