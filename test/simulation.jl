@testset "grid" begin
    max_x = 3e-6
    max_y = 10e-6
    max_t = 5e-12

    nx = 60
    ny = 200

    grid = Grid(max_x, max_y, max_t, nx, ny)

    @test size(grid) == (nx, ny)
    @test size(grid, 1) == nx
    @test size(grid, 2) == ny
    @test axes(grid) == (Base.OneTo(nx), Base.OneTo(ny))
    @test axes(grid, 1) == Base.OneTo(nx)
    @test axes(grid, 2) == Base.OneTo(ny)
    @test boundary(grid) == (max_x, max_y)
end

@testset "Light" begin
    λ = 2.04e-6

    light = Light(λ)

    @test light.λ == λ
    @test light.k == 2π/λ
end

@testset "permittivity" begin
    max_x = 3e-6
    max_y = 10e-6
    max_t = 5e-12

    nx = 60
    ny = 200

    ϵ = 9.

    grid = Grid(max_x, max_y, max_t, nx, ny)
    permittivity = Permittivity(ϵ, grid)

    @test permittivity.ϵ[1] == ϵ
    @test permittivity.ϵx[1] == C * grid.Δt/grid.Δx / ϵ
    @test permittivity.ϵy[1] == C * grid.Δt/grid.Δy / ϵ
    @test all_eq(permittivity.ϵ)
    @test all_eq(permittivity.ϵx)
    @test all_eq(permittivity.ϵy)
end

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

    s = Simulator(nx=nx, ny=ny)
    simulate!(s)

    s = Simulator(grid, light, permittivity, permeability)
    simulate!(s)
end
