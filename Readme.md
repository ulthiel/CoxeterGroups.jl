# CoxeterGroups.jl

[![][docs-dev-img]][docs-dev-url]
[![][action-img]][action-url]
[![][codecov-img]][codecov-url]

A Julia package for computing in Coxeter groups (draft, under development).

The current code is due to T. Schmit (2021) and implements the algorithm from the paper

Casselman, B. (2002). Computation in Coxeter groups. I. Multiplication. *Electron. J. Combin., 9*(1), Research Paper 25, 22.

## Installation

```
julia> using Pkg
julia> Pkg.add(url="https://github.com/ulthiel/CoxeterGroups.jl")
```

## Usage

```julia
julia> using CoxeterGroups
julia> G,(a,b,c) = CoxeterGroup([1 4 2; 4 1 3; 2 3 1], ["a", "b", "c"]);
julia> w = a*c
ca
julia> w*c
a
```

## Todo

See [here](https://github.com/ulthiel/CoxeterGroups.jl/issues/1).


## Development

Checkout the code, change directory in, and start the julia shell.
To run tests, first activate the package, then run the `test` command repeatedly while editing.

    pkg> activate .
    pkg> test

To build the documentation:

    julia --project=docs/  docs/make.jl

You can run a local server to preview the documentation (some of the hyperlinks don't work otherwise):

    python3 -m http.server --directory docs/build/ 8181

and now open <http://localhost:8181/> in the browser.
If the documentation build is complaining about missing dependencies, then run

    $ julia --project=docs/
    julia> ]
    pkg> resolve

<!-- URLS -->


[action-img]: https://github.com/ulthiel/CoxeterGroups.jl/actions/workflows/CI.yml/badge.svg
[action-url]: https://github.com/ulthiel/CoxeterGroups.jl/actions
[codecov-img]: https://codecov.io/gh/ulthiel/CoxeterGroups.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/ulthiel/CoxeterGroups.jl?branch=master
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ulthiel.github.io/CoxeterGroups.jl/dev/
