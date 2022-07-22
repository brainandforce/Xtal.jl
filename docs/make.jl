using Documenter
using Xtal

makedocs(
    sitename = "Xtal",
    format = Documenter.HTML(),
    modules = [Xtal],
    pages = [
        "Home" => "index.md",
        "Getting started" => Any[
            "manual/getting-started.md"
        ],
        "API" => Any[
            "Lattices" => "api/lattices.md"
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
