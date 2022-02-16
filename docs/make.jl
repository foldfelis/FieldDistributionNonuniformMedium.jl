using FieldDistributionNonuniformMedium
using Documenter

DocMeta.setdocmeta!(FieldDistributionNonuniformMedium, :DocTestSetup, :(using FieldDistributionNonuniformMedium); recursive=true)

makedocs(;
    modules=[FieldDistributionNonuniformMedium],
    authors="JingYu Ning <foldfelis@gmail.com> and contributors",
    repo="https://github.com/foldfelis/FieldDistributionNonuniformMedium.jl/blob/{commit}{path}#{line}",
    sitename="FieldDistributionNonuniformMedium.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://foldfelis.github.io/FieldDistributionNonuniformMedium.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Simulation" => "simulation.md",
        "APIs" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/foldfelis/FieldDistributionNonuniformMedium.jl",
    devbranch="master",
)
