using OptimaKit
using Documenter

DocMeta.setdocmeta!(
    OptimaKit,
    :DocTestSetup,
    :(using OptimaKit; import OptimaKit: solve);
    recursive = true,
)

makedocs(;
    clean    = false,
    modules  = [OptimaKit],
    authors  = "Jean-François Barthélémy",
    sitename = "OptimaKit.jl",
    remotes  = Dict(".." => Documenter.Remotes.GitHub("ChemistryTools", "OptimaKit.jl")),
    format   = Documenter.HTML(;
        canonical  = "https://ChemistryTools.github.io/OptimaKit.jl",
        edit_link  = "main",
        prettyurls = (get(ENV, "CI", nothing) == "true"),
        collapselevel = 1,
    ),
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "Theory"          => "theory.md",
        "Examples"        => [
            "Basic Usage"     => "examples/basic_usage.md",
            "Warm Start"      => "examples/warm_start.md",
            "Sensitivity"     => "examples/sensitivity.md",
            "SciML Interface" => "examples/sciml_interface.md",
        ],
        "API Reference"   => "api.md",
    ],
    warnonly = [:missing_docs, :docs_block],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(; repo = "github.com/ChemistryTools/OptimaKit.jl", devbranch = "main")
end
