using BattPhase
using Documenter

DocMeta.setdocmeta!(BattPhase, :DocTestSetup, :(using BattPhase); recursive=true)

makedocs(;
    modules=[BattPhase],
    authors="Brady Planden",
    repo="https://github.com/bradyplanden/BattPhase.jl/blob/{commit}{path}#{line}",
    sitename="BattPhase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bradyplanden.github.io/BattPhase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bradyplanden/BattPhase.jl",
    devbranch="main",
)
