using Documenter

println(LOAD_PATH)
using SARM

mathengine = Documenter.Writers.HTMLWriter.KaTeX()
#mathengine = Documenter.Writers.HTMLWriter.MathJax2()

mathengine = MathJax3(Dict(
    :loader => Dict("load" => ["[tex]/physics"]),
    :tex => Dict(
        "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
        "tags" => "ams",
        "packages" => ["base", "ams", "autoload", "physics"],
    ),
))

format = Documenter.Writers.HTMLWriter.HTML(mathengine = mathengine, prettyurls = false)
#Documenter.Writers.LaTeXWriter.Latex()

#    makedocs(sitename = "Opt40", authors = "Rolf Wester")
#        modules = [Opt40],
#        clean = false,
#        format = Documenter.Writers.HTMLWriter(),
#        #analytics = "",
#        #linkcheck = !("skiplinks" in ARGS),
#        pages = Any[ # Compat: `Any` for 0.4 compat
#            "Home" => "index.md",
#            ]
#        # Use clean URLs, unless built as a "local" build
#        #html_prettyurls = !("local" in ARGS),
#
#    makedocs(
#        root    = "<current-directory>",
#        source  = "src",
#        build   = "build",
#        clean   = true,
#        doctest = true,
#        modules = Module[],
#        repo    = "",
#        highlightsig = true,
#        sitename = "",
#        expandfirst = [],
#        draft = false,


makedocs(sitename = "SimpleRadTrans", authors = "Rolf Wester", modules = [SimpleRadTrans], format=format)

#julia --color=yes make.jl
