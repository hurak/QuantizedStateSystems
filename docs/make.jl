using QuantizedStateSystems
using Documenter

DocMeta.setdocmeta!(QuantizedStateSystems, :DocTestSetup, :(using QuantizedStateSystems); recursive=true)

makedocs(;
    modules=[QuantizedStateSystems],
    authors="Zdeněk Hurák <hurak@fel.cvut.cz>",
    repo="https://github.com/hurak/QuantizedStateSystems.jl/blob/{commit}{path}#{line}",
    sitename="QuantizedStateSystems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hurak.github.io/QuantizedStateSystems.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hurak/QuantizedStateSystems.jl",
    devbranch="main",
)
