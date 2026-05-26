using Test
using DataFrames
using DoMarketsReduceCostsReplication

@testset "Helper functions" begin
    toy = DataFrame(plant_num2 = [1, 1, 2, 3, 3, 3])

    @test count_plant_epochs(toy) == 3
    @test endswith(raw_data_dir(), joinpath("data", "raw", "openicpsr", "116286-V1"))
end

@testset "Input data preparation" begin
    toy = DataFrame(
        plant_num2 = [1, 1, 1, 2, 2],
        yr_data = [1981, 1982, 1984, 1981, 1982],
        ioudum = [1, 1, 1, 0, 0],
        hi_nug = [1, 1, 1, 1, 1],
    )

    prepared = prepare_input_data(toy)

    @test prepared.diff == [false, true, false, false, true]
    @test prepared.iounug == [0, 0, 0, 0, 0]
end

@testset "Prais transform" begin
    x = [10.0, 12.0, 20.0]
    diff = [false, true, false]
    rho = 0.5

    transformed = prais_transform(x, diff, rho)

    @test transformed[1] ≈ sqrt(1 - rho^2) * 10.0
    @test transformed[2] ≈ 12.0 - rho * 10.0
    @test transformed[3] ≈ sqrt(1 - rho^2) * 20.0
end
@testset "Rho estimation" begin
    errors = [2.0, 1.0, 3.0, 1.5]
    diff = [false, true, false, true]

    rho = estimate_rho(errors, diff)

    expected = (2.0 * 1.0 + 3.0 * 1.5) / (2.0^2 + 3.0^2)
    @test rho ≈ expected
end
@testset "Add Prais columns" begin
    toy = DataFrame(
        diff = [false, true],
        x = [10.0, 12.0],
    )

    add_prais_columns!(toy, [:x], 0.5)

    @test "prais_x" in names(toy)
    @test toy.prais_x[1] ≈ sqrt(1 - 0.5^2) * 10.0
    @test toy.prais_x[2] ≈ 12.0 - 0.5 * 10.0
end