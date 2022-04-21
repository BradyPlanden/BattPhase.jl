using LiMetalPhaseFields
using Documenter

DocMeta.setdocmeta!(LiMetalPhaseFields, :DocTestSetup, :(using LiMetalPhaseFields); recursive=true)

makedocs(;
    modules=[LiMetalPhaseFields],
    authors="Brady Planden",
    repo="https://github.com/bradyplanden/LiMetalPhaseFields.jl/blob/{commit}{path}#{line}",
    sitename="LiMetalPhaseFields.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bradyplanden.github.io/LiMetalPhaseFields.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bradyplanden/LiMetalPhaseFields.jl",
    devbranch="main",
)
