using FieldDistributionNonuniformMedium
using Test

const C = 299792458

all_eq(x::AbstractArray) = all(xᵢ -> xᵢ == x[1], x)

@testset "FieldDistributionNonuniformMedium.jl" begin
    include("simulation.jl")
    include("plot.jl")
end
