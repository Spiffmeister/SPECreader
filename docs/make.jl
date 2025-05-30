
push!(LOAD_PATH,"../src/")

using Documenter, Literate
using Pkg; Pkg.activate(".."); using SPECreader


makedocs(sitename="SPECreader",
    pages = [
        "Home" => "index.md",
        "Examples" => [
        ]
    ],
    modules=[SPECreader],
    format=Documenter.HTML(prettyurls=false),
    warnonly = Documenter.except(:linkcheck,:footnote)
    )

deploydocs(
    repo = "github.com/Spiffmeister/SPECreader.jl.git",
)