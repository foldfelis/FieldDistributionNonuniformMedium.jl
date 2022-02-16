var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FieldDistributionNonuniformMedium","category":"page"},{"location":"#FieldDistributionNonuniformMedium","page":"Home","title":"FieldDistributionNonuniformMedium","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FieldDistributionNonuniformMedium.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/foldfelis/FieldDistributionNonuniformMedium.jl","category":"page"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using FieldDistributionNonuniformMedium\n\njulia> s = Simulator();\n\njulia> simulate!(s);\n\njulia> plot_e_field(s)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#APIs","page":"Home","title":"APIs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [FieldDistributionNonuniformMedium]","category":"page"},{"location":"#FieldDistributionNonuniformMedium.Grid-NTuple{5, Any}","page":"Home","title":"FieldDistributionNonuniformMedium.Grid","text":"Grid(max_x, max_y, max_t, nx, ny)\n\nConstruct a discretization grid for computational domain.\n\nArguments\n\nmax_x: Linear sizes of first computational domain in meters.\nmax_y: Linear sizes of second computational domain in meters.\nmax_t: Maximum computational time in seconds.\nnx: Number of discretization grid for first computational domain.\nny: Number of discretization grid for second computational domain.\n\nExample\n\njulia> Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\n\n\n\n\n","category":"method"},{"location":"#FieldDistributionNonuniformMedium.Light-Tuple{Any}","page":"Home","title":"FieldDistributionNonuniformMedium.Light","text":"Light(λ)\n\nConstruct a light to propagate.\n\nArguments\n\nλ: Wavelength of the light in meters.\n\nExample\n\njulia> Light(2.04e-6);\n\n\n\n\n\n","category":"method"},{"location":"#FieldDistributionNonuniformMedium.Permeability-Tuple{Real, Grid}","page":"Home","title":"FieldDistributionNonuniformMedium.Permeability","text":"Permeability(μ_const::Real, grid::Grid)\n\nConstruct a uniform permeability for the medium.\n\nArguments\n\nμ_const: Permeability of the medium in N/A².\ngrid: Discretization grid.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> Permeability(1., grid);\n\n\n\n\n\n","category":"method"},{"location":"#FieldDistributionNonuniformMedium.Permittivity-Tuple{Real, Grid}","page":"Home","title":"FieldDistributionNonuniformMedium.Permittivity","text":"Permittivity(ϵ_const::Real, grid::Grid)\n\nConstruct a uniform permittivity for the medium.\n\nArguments\n\nϵ_const: Permittivity of the medium in F/m.\ngrid: Discretization grid.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> Permittivity(9., grid);\n\n\n\n\n\n","category":"method"},{"location":"#FieldDistributionNonuniformMedium.Simulator-Tuple{Grid, Light, Permittivity, Permeability}","page":"Home","title":"FieldDistributionNonuniformMedium.Simulator","text":"Simulator(grid::Grid, light::Light, permittivity::Permittivity, permeability::Permeability)\n\nArguments\n\ngrid: Discretization grid for computational domain.\nlight: Light to propagate.\npermittivity: Permittivity for the medium.\npermeability: Permeability for the medium.\n\nExample\n\n```jldoctest julia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> light = Light(2.04e-6);\n\njulia> permittivity = Permittivity(9., grid);\n\njulia> permeability = Permeability(1., grid);\n\njulia> Simulator(grid, light, permittivity, permeability);\n\n\n\n\n\n","category":"method"},{"location":"#FieldDistributionNonuniformMedium.implant!-Tuple{Permittivity, Real, AbstractVector, AbstractVector, AbstractVector, Grid}","page":"Home","title":"FieldDistributionNonuniformMedium.implant!","text":"implant!(\n    permittivity::Permittivity, ϵ_const::Real,\n    xs::AbstractVector, ys::AbstractVector, rs::AbstractVector,\n    grid::Grid\n)\n\nImplant some bubble defect into the medium.\n\nArguments\n\npermittivity: Permittivity object.\nϵ_const: Permittivity of the defect in F/m.\nxs: Position of defect of first computational domain in meters.\nys: Position of defect of second computational domain in meters.\nrs: Radius of defect in meters.\ngrid: Discretization grid.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> permittivity = Permittivity(9., grid);\n\njulia> ϵ_defect = 1.;\n\njulia> xs_defect = [0, 1e-6, -1e-6];\n\njulia> ys_defect = [1e-6, 2e-6, 3e-6];\n\njulia> rs_defect = [0.5e-6, 0.1e-6, 0.2e-6];\n\njulia> implant!(\n           permittivity, ϵ_defect,\n           xs_defect, ys_defect, rs_defect,\n           grid\n       );\n\n\n\n\n\n","category":"method"},{"location":"#FieldDistributionNonuniformMedium.simulate!-Tuple{Simulator}","page":"Home","title":"FieldDistributionNonuniformMedium.simulate!","text":"simulate!(s::Simulator)\n\nRun simulation from current t to max_t\n\nArguments\n\ns: Simulator.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> light = Light(2.04e-6);\n\njulia> permittivity = Permittivity(9., grid);\n\njulia> permeability = Permeability(1., grid);\n\njulia> s = Simulator(grid, light, permittivity, permeability);\n\njulia> simulate!(s);\n\n\n\n\n\n","category":"method"}]
}
