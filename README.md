# FieldDistributionNonuniformMedium

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://foldfelis.github.io/FieldDistributionNonuniformMedium.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://foldfelis.github.io/FieldDistributionNonuniformMedium.jl/dev)
[![Build Status](https://github.com/foldfelis/FieldDistributionNonuniformMedium.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/foldfelis/FieldDistributionNonuniformMedium.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/foldfelis/FieldDistributionNonuniformMedium.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/foldfelis/FieldDistributionNonuniformMedium.jl)

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add https://github.com/foldfelis/FieldDistributionNonuniformMedium.jl
```

## Quick start

```julia
julia> using FieldDistributionNonuniformMedium

julia> s = Simulator();

julia> simulate!(s);

julia> plot_e_field(s)
```

![](docs/src/assets/demo.png)
