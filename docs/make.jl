using Documenter
using MixedModelsMakie

makedocs(;
         sitename="MixedModelsMakie",
         doctest=true,
         checkdocs=:exports,
         format=Documenter.HTML(; 
                                size_threshold_ignore=["api.md"]),
         pages=["index.md",
                "api.md"])

deploydocs(; repo="github.com/JuliaMixedModels/MixedModelsMakie.jl.git", devbranch="main",
           push_preview=true)
