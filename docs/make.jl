import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))


using Documenter, Chion

makedocs(sitename="My Documentation")
