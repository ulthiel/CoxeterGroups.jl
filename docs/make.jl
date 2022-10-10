using Documenter
using CoxeterGroups
import CoxeterGroups: CoxeterSystem

DocMeta.setdocmeta!(CoxeterGroups, :DocTestSetup, :(using CoxeterGroups); recursive=true)

makedocs(;
    modules=[CoxeterGroups],
    sitename="CoxeterGroups.jl"
)

deploydocs(;
    repo="github.com/ulthiel/CoxeterGroups.jl",
)