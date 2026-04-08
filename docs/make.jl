using OptimaJL
using Documenter

DocMeta.setdocmeta!(
    OptimaJL,
    :DocTestSetup,
    :(using OptimaJL; import OptimaJL: solve);
    recursive = true,
)

makedocs(;
    clean    = false,
    modules  = [OptimaJL],
    authors  = "Jean-François Barthélémy",
    sitename = "OptimaJL.jl",
    remotes  = Dict(".." => Documenter.Remotes.GitHub("jfbarthelemy", "OptimaJL.jl")),
    format   = Documenter.HTML(;
        canonical  = "https://jfbarthelemy.github.io/OptimaJL.jl",
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
    deploydocs(; repo = "github.com/jfbarthelemy/OptimaJL.jl", devbranch = "main")
end
