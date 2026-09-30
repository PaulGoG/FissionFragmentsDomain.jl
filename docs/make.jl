include(joinpath(@__DIR__, "activate.jl"))

using FissionFragmentsDomain
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :numeric)

DocMeta.setdocmeta!(
    FissionFragmentsDomain,
    :DocTestSetup,
    :(using FissionFragmentsDomain);
    recursive = true,
)

makedocs(;
    modules = [FissionFragmentsDomain],
    authors = "Paul-Adrian Gogîță",
    sitename = "FissionFragmentsDomain.jl",
    repo = Documenter.Remotes.GitHub("PaulGoG", "FissionFragmentsDomain.jl"),
    format = Documenter.HTML(;
        canonical = "https://PaulGoG.github.io/FissionFragmentsDomain.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Systems, tables and masses" => "system.md",
        "Charge distribution and domain" => "fragmentation.md",
        "Level density and energetics" => "energetics.md",
        "Yields and averages" => "yields.md",
        "Curves and run records" => "curves.md",
        "References" => "references.md",
    ],
    plugins = [bib],
)

deploydocs(; repo = "github.com/PaulGoG/FissionFragmentsDomain.jl", devbranch = "main")
