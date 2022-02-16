export
    Grid,
    build,
    boundary,

    Light,

    Permittivity,
    implant!,

    Permeability,

    Simulator,
    next!,
    simulate!

const C = 299792458

struct Grid{T<:Real}
    nx::Int
    ny::Int
    nt::Int

    Δx::T
    Δy::T
    Δt::T

    max_x::T
    max_y::T
    max_t::T
end

function Grid(max_x, max_y, max_t, nx, ny)
    Δx = max_x / nx
    Δy = max_y / ny

    Δt = 1 / C / √(1/Δx^2 + 1/Δy^2)
    nt = round(Int, max_t/Δt)

    return Grid(nx, ny, nt, Δx, Δy, Δt, max_x, max_y, max_t)
end

function build(grid::Grid)
    return cat(
        repeat(axes(grid, 1), 1, size(grid, 2), 1),
        repeat(axes(grid, 2)', size(grid, 1), 1, 1),
        dims=3
    )
end

Base.size(grid::Grid) = (grid.nx, grid.ny)
Base.size(grid::Grid, d) = d::Integer <= 2 ? size(grid)[d] : 1
function Base.axes(grid::Grid)
    axes_x = grid.Δx * (Base.OneTo(grid.nx) .- grid.nx/2)
    axes_y = grid.Δy * (Base.OneTo(grid.ny) .- grid.Δy/2)

    return (axes_x, axes_y)
end
Base.axes(grid::Grid, d) = d::Integer <= 2 ? axes(grid)[d] : 1
boundary(grid::Grid) = (grid.max_x, grid.max_y)
boundary(grid::Grid, d) = d::Integer <= 2 ? boundary(grid)[d] : 1

struct Light{T<:Real}
    λ::T
    k::T
end

function Light(λ)
    return Light(λ, 2π/λ)
end

struct Permittivity{T<:AbstractMatrix}
    ϵ::T
    ϵx::T
    ϵy::T
end

function Permittivity(ϵ_const::Real, grid::Grid)
    ϵ = ϵ_const * ones(size(grid))
    ϵx = C * grid.Δt/grid.Δx ./ ϵ
    ϵy = C * grid.Δt/grid.Δy ./ ϵ

    return Permittivity(ϵ, ϵx, ϵy)
end

function implant!(permittivity::Permittivity, ϵ_const::Real, xs::AbstractVector, ys::AbstractVector, rs::AbstractVector, grid::Grid)
    length(xs) == length(ys) == length(rs) || throw(DimensionMismatch("xs, ys, rs must have same length"))

    in_circle(i, j) = true in [
        √(
            (axes(grid, 1)[i] - xs[c])^2 +
            (axes(grid, 2)[j] - ys[c])^2
        ) < rs[c]
        for c in 1:length(xs)
    ]

    for i in 1:size(grid, 1), j in 1:size(grid, 2)
        in_circle(i, j) && (permittivity.ϵ[i, j] = ϵ_const)
    end

    permittivity.ϵx .= C * grid.Δt/grid.Δx ./ permittivity.ϵ
    permittivity.ϵy .= C * grid.Δt/grid.Δy ./ permittivity.ϵ

    return permittivity
end

# function Base.rand(Tϵ::Type{Permittivity}, n::Integer, r::Real, grid::Grid)
#     xs = grid.max_x .* rand(n)
#     ys = grid.max_y .* rand(n)
#     rs = r .* rand(n)

#     return implant(Tϵ, xs, ys, rs, grid)
# end

struct Permeability{T<:AbstractMatrix}
    μ::T
    μx::T
    μy::T
end

function Permeability(μ_const::Real, grid::Grid)
    μ = μ_const * ones(size(grid))
    μx = C * grid.Δt/grid.Δx ./ μ
    μy = C * grid.Δt/grid.Δy ./ μ

    return Permeability(μ, μx, μy)
end

mutable struct Simulator{T<:AbstractArray}
    grid::Grid
    light::Light
    permittivity::Permittivity
    permeability::Permeability

    ez::T
    hx::T
    hy::T

    t::Int
end

function get_default_e_field(light::Light, grid::Grid; t=1)
    Δt = grid.Δt
    k = light.k

    return 0.1 * exp.(
        -axes(grid, 1)[2:end].^2 ./
        (boundary(grid, 1)/4)^2
    ) * sin(k * C*Δt*t)
end

function Simulator(grid::Grid, light::Light, permittivity::Permittivity, permeability::Permeability)
    ez = zeros(Float64, size(grid))
    ez[2:end, 1] .= get_default_e_field(light, grid)

    return Simulator(
        grid,
        light,
        permittivity,
        permeability,

        ez,
        zeros(Float64, size(grid)),
        zeros(Float64, size(grid)),

        1
    )
end

function Simulator(;
    max_x=3e-6, max_y=10e-6, max_t=5e-12,
    nx=300, ny=1000,
    λ=2.04e-6,
    ϵ = 9., μ = 1.,
    ndefect=rand(1:5), r=0.45e-6,
)
    grid = Grid(max_x, max_y, max_t, nx, ny)
    light = Light(λ)
    permittivity = Permittivity(ϵ, grid)
    permeability = Permeability(μ, grid)

    return Simulator(grid, light, permittivity, permeability)
end

function next!(s::Simulator)
    s.t += 1

    ϵx, ϵy = s.permittivity.ϵx, s.permittivity.ϵy
    μx, μy = s.permeability.μx, s.permeability.μy

    s.ez[2:end, 1] .+= get_default_e_field(s.light, s.grid, t=s.t)

    s.hx[2:end-1, 2:end-1] .+= -μx[2:end-1, 2:end-1].*(s.ez[2:end-1, 2:end-1] - s.ez[2:end-1, 1:end-2])
    s.hy[2:end-1, 2:end-1] .+= +μy[2:end-1, 2:end-1].*(s.ez[2:end-1, 2:end-1] - s.ez[1:end-2, 2:end-1])

    s.ez[2:end-1, 2:end-1] .+=
        ϵx[2:end-1, 2:end-1].*(s.hy[3:end, 2:end-1] - s.hy[2:end-1, 2:end-1]) -
        ϵy[2:end-1, 2:end-1].*(s.hx[2:end-1, 3:end] - s.hx[2:end-1, 2:end-1])

    return s
end

function simulate!(s::Simulator)
    nt = s.grid.nt

    p = Progress(nt)
    for _ in 2:nt
        next!(s)
        ProgressMeter.next!(p)
    end

    return s
end
