@testset "grid" begin
    max_x = 3e-6
    max_y = 10e-6
    max_t = 5e-12

    nx = 60
    ny = 200

    Δx = max_x / nx
    Δy = max_y / ny

    Δt = 1 / C / √(1/Δx^2 + 1/Δy^2)
    nt = round(Int, max_t/Δt)

    axes_x = Δx * (Base.OneTo(nx) .- nx/2)
    axes_y = Δy * (Base.OneTo(ny) .- Δy/2)

    grid = Grid(max_x, max_y, max_t, nx, ny)

    @test grid.Δx == Δx
    @test grid.Δy == Δy
    @test grid.Δt == Δt
    @test grid.nt == nt

    @test size(grid) == (nx, ny)
    @test size(grid, 1) == nx
    @test size(grid, 2) == ny
    @test axes(grid) == (axes_x, axes_y)
    @test axes(grid, 1) == axes_x
    @test axes(grid, 2) == axes_y
    @test boundary(grid) == (max_x, max_y)
    @test boundary(grid, 1) == max_x
    @test boundary(grid, 2) == max_y
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

@testset "permeability" begin
    max_x = 3e-6
    max_y = 10e-6
    max_t = 5e-12

    nx = 60
    ny = 200

    μ = 1.

    grid = Grid(max_x, max_y, max_t, nx, ny)
    permeability = Permeability(μ, grid)

    @test permeability.μ[1] == μ
    @test permeability.μx[1] == C * grid.Δt/grid.Δx / μ
    @test permeability.μy[1] == C * grid.Δt/grid.Δy / μ
    @test all_eq(permeability.μ)
    @test all_eq(permeability.μx)
    @test all_eq(permeability.μy)
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
