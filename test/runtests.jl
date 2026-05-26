using Test
using DataFrames
using DoMarketsReduceCostsReplication

@testset "Helper functions" begin
    toy = DataFrame(plant_num2 = [1, 1, 2, 3, 3, 3])

    @test count_plant_epochs(toy) == 3
    @test endswith(raw_data_dir(), joinpath("data", "raw", "openicpsr", "116286-V1"))
end