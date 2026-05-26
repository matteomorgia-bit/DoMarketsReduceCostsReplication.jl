using Test
using DataFrames
using StatFiles

const RAW_DIR = joinpath(@__DIR__, "..", "data", "raw", "openicpsr", "116286-V1")

function load_stata_file(filename)
    path = joinpath(RAW_DIR, filename)
    return DataFrame(load(path))
end

@testset "FRW estimation samples" begin
    enf = load_stata_file("frw1extract_enf.dta")
    fuel = load_stata_file("frw1extract_f.dta")

    @test nrow(enf) == 10079
    @test length(unique(enf.plant_num2)) == 769

    @test nrow(fuel) == 10002
    @test length(unique(fuel.plant_num2)) == 768
end 