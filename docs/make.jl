using Documenter
using MixedModelsMakie

makedocs(;
         sitename="MixedModelsMakie",
         doctest=true,
         checkdocs=:exports,
         pages=["index.md",
                "api.md"])

deploydocs(; repo="github.com/palday/MixedModelsMakie.jl.git", devbranch="main",
           push_preview=true)
