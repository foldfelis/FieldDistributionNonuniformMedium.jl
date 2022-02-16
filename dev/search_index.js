var documenterSearchIndex = {"docs":
[{"location":"api/#Index","page":"APIs","title":"Index","text":"","category":"section"},{"location":"api/","page":"APIs","title":"APIs","text":"","category":"page"},{"location":"api/#APIs","page":"APIs","title":"APIs","text":"","category":"section"},{"location":"api/","page":"APIs","title":"APIs","text":"Modules = [FieldDistributionNonuniformMedium]","category":"page"},{"location":"api/#FieldDistributionNonuniformMedium.Grid-NTuple{5, Any}","page":"APIs","title":"FieldDistributionNonuniformMedium.Grid","text":"Grid(max_x, max_y, max_t, nx, ny)\n\nConstruct a discretization grid for computational domain.\n\nArguments\n\nmax_x: Linear sizes of first computational domain in meters.\nmax_y: Linear sizes of second computational domain in meters.\nmax_t: Maximum computational time in seconds.\nnx: Number of discretization grid for first computational domain.\nny: Number of discretization grid for second computational domain.\n\nExample\n\njulia> Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\n\n\n\n\n","category":"method"},{"location":"api/#FieldDistributionNonuniformMedium.Light-Tuple{Any}","page":"APIs","title":"FieldDistributionNonuniformMedium.Light","text":"Light(λ)\n\nConstruct a light to propagate.\n\nArguments\n\nλ: Wavelength of the light in meters.\n\nExample\n\njulia> Light(2.04e-6);\n\n\n\n\n\n","category":"method"},{"location":"api/#FieldDistributionNonuniformMedium.Permeability-Tuple{Real, Grid}","page":"APIs","title":"FieldDistributionNonuniformMedium.Permeability","text":"Permeability(μ_const::Real, grid::Grid)\n\nConstruct a uniform permeability for the medium.\n\nArguments\n\nμ_const: Permeability of the medium in N/A².\ngrid: Discretization grid.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> Permeability(1., grid);\n\n\n\n\n\n","category":"method"},{"location":"api/#FieldDistributionNonuniformMedium.Permittivity-Tuple{Real, Grid}","page":"APIs","title":"FieldDistributionNonuniformMedium.Permittivity","text":"Permittivity(ϵ_const::Real, grid::Grid)\n\nConstruct a uniform permittivity for the medium.\n\nArguments\n\nϵ_const: Permittivity of the medium in F/m.\ngrid: Discretization grid.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> Permittivity(9., grid);\n\n\n\n\n\n","category":"method"},{"location":"api/#FieldDistributionNonuniformMedium.Simulator-Tuple{Grid, Light, Permittivity, Permeability}","page":"APIs","title":"FieldDistributionNonuniformMedium.Simulator","text":"Simulator(grid::Grid, light::Light, permittivity::Permittivity, permeability::Permeability)\n\nArguments\n\ngrid: Discretization grid for computational domain.\nlight: Light to propagate.\npermittivity: Permittivity for the medium.\npermeability: Permeability for the medium.\n\nExample\n\n```jldoctest julia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> light = Light(2.04e-6);\n\njulia> permittivity = Permittivity(9., grid);\n\njulia> permeability = Permeability(1., grid);\n\njulia> Simulator(grid, light, permittivity, permeability);\n\n\n\n\n\n","category":"method"},{"location":"api/#FieldDistributionNonuniformMedium.implant!-Tuple{Permittivity, Real, AbstractVector, AbstractVector, AbstractVector, Grid}","page":"APIs","title":"FieldDistributionNonuniformMedium.implant!","text":"implant!(\n    permittivity::Permittivity, ϵ_const::Real,\n    xs::AbstractVector, ys::AbstractVector, rs::AbstractVector,\n    grid::Grid\n)\n\nImplant some bubble defect into the medium.\n\nArguments\n\npermittivity: Permittivity object.\nϵ_const: Permittivity of the defect in F/m.\nxs: Position of defect of first computational domain in meters.\nys: Position of defect of second computational domain in meters.\nrs: Radius of defect in meters.\ngrid: Discretization grid.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> permittivity = Permittivity(9., grid);\n\njulia> ϵ_defect = 1.;\n\njulia> xs_defect = [0, 1e-6, -1e-6];\n\njulia> ys_defect = [1e-6, 2e-6, 3e-6];\n\njulia> rs_defect = [0.5e-6, 0.1e-6, 0.2e-6];\n\njulia> implant!(\n           permittivity, ϵ_defect,\n           xs_defect, ys_defect, rs_defect,\n           grid\n       );\n\n\n\n\n\n","category":"method"},{"location":"api/#FieldDistributionNonuniformMedium.simulate!-Tuple{Simulator}","page":"APIs","title":"FieldDistributionNonuniformMedium.simulate!","text":"simulate!(s::Simulator)\n\nRun simulation from current t to max_t\n\nArguments\n\ns: Simulator.\n\nExample\n\njulia> grid = Grid(3e-6, 10e-6, 1e-12, 60, 200);\n\njulia> light = Light(2.04e-6);\n\njulia> permittivity = Permittivity(9., grid);\n\njulia> permeability = Permeability(1., grid);\n\njulia> s = Simulator(grid, light, permittivity, permeability);\n\njulia> simulate!(s);\n\n\n\n\n\n","category":"method"},{"location":"simulation/#Simulation","page":"Simulation","title":"Simulation","text":"","category":"section"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"There are two way to construct a Simulator with uniform permittivity.","category":"page"},{"location":"simulation/#Construct-simulator-by-parameter","page":"Simulation","title":"Construct simulator by parameter","text":"","category":"section"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"s = Simulator(\n    max_x=3e-6, max_y=10e-6, max_t=0.1e-12, # calculation boundaries\n    nx=300, ny=1000, # discretization\n    λ=2.04e-6, # wavelength of light\n    ϵ = 9., μ = 1. # permittivity and permeability\n)","category":"page"},{"location":"simulation/#Construct-simulator-by-components","page":"Simulation","title":"Construct simulator by components","text":"","category":"section"},{"location":"simulation/#Declare-constants","page":"Simulation","title":"Declare constants","text":"","category":"section"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"# calculation boundaries\nmax_x = 3e-6\nmax_y = 10e-6\nmax_t = 0.1e-12\n\n# discretization\nnx = 300\nny = 1000\n\n# wavelength of light\nλ = 2.04e-6\n\n# permittivity and permeability\nϵ = 9.\nμ = 1.","category":"page"},{"location":"simulation/#Construct-components","page":"Simulation","title":"Construct components","text":"","category":"section"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"grid = Grid(max_x, max_y, max_t, nx, ny)\nlight = Light(λ)\npermittivity = Permittivity(ϵ, grid)\npermeability = Permeability(μ, grid)","category":"page"},{"location":"simulation/#Construct-simulator","page":"Simulation","title":"Construct simulator","text":"","category":"section"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"s = Simulator(grid, light, permittivity, permeability)","category":"page"},{"location":"simulation/#Implant-defect-to-modify-permittivity","page":"Simulation","title":"Implant defect to modify permittivity","text":"","category":"section"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"# ##########\n# # const. #\n# ##########\n\n# calculation boundaries\nmax_x = 3e-6\nmax_y = 10e-6\nmax_t = 0.1e-12\n\n# discretization\nnx = 300\nny = 1000\n\n# wavelength of light\nλ = 2.04e-6\n\n# permittivity and permeability\nϵ = 9.\nμ = 1.\n\n# defect\nϵ_defect = 1.\nxs_defect = [0, 1e-6, -1e-6]\nys_defect = [1e-6, 2e-6, 3e-6]\nrs_defect = [0.5e-6, 0.1e-6, 0.2e-6]\n\n# ##############\n# # components #\n# ##############\n\ngrid = Grid(max_x, max_y, max_t, nx, ny)\nlight = Light(λ)\npermittivity = Permittivity(ϵ, grid)\npermeability = Permeability(μ, grid)\n\nimplant!(permittivity, ϵ_defect, xs_defect, ys_defect, rs_defect, grid)\n\n# #############\n# # simulator #\n# #############\ns = Simulator(grid, light, permittivity, permeability)","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"Or","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"ϵ_defect = 1.\nxs_defect = [0, 1e-6, -1e-6]\nys_defect = [1e-6, 2e-6, 3e-6]\nrs_defect = [0.5e-6, 0.1e-6, 0.2e-6]\n\ns = Simulator(\n    max_x=3e-6, max_y=10e-6, max_t=0.1e-12, # calculation boundaries\n    nx=300, ny=1000, # discretization\n    λ=2.04e-6, # wavelength of light\n    ϵ = 9., μ = 1. # permittivity and permeability\n)\n\nimplant!(s.permittivity, ϵ_defect, xs_defect, ys_defect, rs_defect, s.grid)","category":"page"},{"location":"simulation/#Run-simulation","page":"Simulation","title":"Run simulation","text":"","category":"section"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"simulate!(s)","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"To see permittivity:","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"plot_ϵ(s)","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"(Image: )","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"To sea the result:","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"plot_e_field(s)","category":"page"},{"location":"simulation/","page":"Simulation","title":"Simulation","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FieldDistributionNonuniformMedium","category":"page"},{"location":"#FieldDistributionNonuniformMedium","page":"Home","title":"FieldDistributionNonuniformMedium","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FieldDistributionNonuniformMedium.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/foldfelis/FieldDistributionNonuniformMedium.jl","category":"page"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using FieldDistributionNonuniformMedium\n\njulia> s = Simulator();\n\njulia> simulate!(s);\n\njulia> plot_e_field(s)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"}]
}
