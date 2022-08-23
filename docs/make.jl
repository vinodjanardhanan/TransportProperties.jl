using TransportProperties
using Documenter

DocMeta.setdocmeta!(TransportProperties, :DocTestSetup, :(using TransportProperties); recursive=true)

makedocs(;
    modules=[TransportProperties],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/TransportProperties.jl/blob/{commit}{path}#{line}",
    sitename="TransportProperties.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/TransportProperties.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/TransportProperties.jl",
    devbranch="main",
)
