@testset "simulation" begin
    # ##########
    # # const. #
    # ##########

    max_x = 3e-6
    max_y = 10e-6
    max_t = 5e-12

    nx = 60
    ny = 200

    λ = 2.04e-6

    ϵ = 9.
    μ = 1.

    ndefect = 5
    r = 0.45e-6

    # ##############
    # # components #
    # ##############

    grid = Grid(max_x, max_y, max_t, nx, ny)
    light = Light(λ)
    permittivity = Permittivity(ϵ, grid)
    permeability = Permeability(μ, grid)

    # #############
    # # simulator #
    # #############

    s = Simulator(grid, light, permittivity, permeability)
    simulate!(s)
end
