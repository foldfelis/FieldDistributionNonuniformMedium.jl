using VisualRegressionTests
using Plots

@testset "utils" begin
    gr()
    
    # ##########
    # # const. #
    # ##########

    max_x = 3e-6
    max_y = 10e-6
    max_t = 0.1e-12

    nx = 60
    ny = 200

    λ = 2.04e-6

    ϵ = 9.
    μ = 1.

    ϵ_defect = 1.
    xs_defect = [0, 1e-6, -1e-6]
    ys_defect = [1e-6, 2e-6, 3e-6]
    rs_defect = [0.5e-6, 0.1e-6, 0.2e-6]

    # ##############
    # # components #
    # ##############

    grid = Grid(max_x, max_y, max_t, nx, ny)
    light = Light(λ)
    permittivity = Permittivity(ϵ, grid)
    permeability = Permeability(μ, grid)

    implant!(permittivity, ϵ_defect, xs_defect, ys_defect, rs_defect, grid)

    # #############
    # # simulator #
    # #############
    s = Simulator(grid, light, permittivity, permeability)
    simulate!(s)

    @plottest begin
        plot_ϵ(s)
    end "assets/permittivity.png"

    @plottest begin
        plot_e_field(s)
    end "assets/e_field.png"
end
