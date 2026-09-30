using FissionFragmentsDomain
using Documenter

DocMeta.setdocmeta!(FissionFragmentsDomain, :DocTestSetup, :(using FissionFragmentsDomain); recursive=true)

makedocs(;
    modules=[FissionFragmentsDomain],
    authors="Paul-Adrian Gogîță",
    sitename="FissionFragmentsDomain.jl",
    format=Documenter.HTML(;
        canonical="https://PaulGoG.github.io/FissionFragmentsDomain.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PaulGoG/FissionFragmentsDomain.jl",
    devbranch="main",
)
