module DoMarketsReduceCostsReplication

export raw_data_dir

"""
    raw_data_dir()

Return the path to the openICPSR raw data folder used by the replication scripts.
"""
function raw_data_dir()
    return joinpath(@__DIR__, "..", "data", "raw", "openicpsr", "116286-V1")
end

end