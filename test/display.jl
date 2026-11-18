@testset "Display" begin
    @test sprint(show, Chion.ThreadsBackend()) == "threads"
    @test sprint(show, Chion.GPUBackend()) == "gpu"

    grid = SnowpackGrid(2)
    forcing = SnowpackForcing(dt_days = [1.0, 1.0], ncol = 2, air_temperature_c = [-5.0, -4.0],
        snowfall_mm_day = 1.0, rainfall_mm_day = 0.0, shortwave_down = [100.0, 110.0], wind_speed = 3.0,
        latitude_deg = 70.0, time_values = [DateTime(2001, 1, 1, 12), DateTime(2001, 1, 2, 12)])
    sim = Simulation(BESSIModel(grid); forcing)
    text = sprint(show, MIME"text/plain"(), sim)
    @test occursin("Simulation", text)
    @test !isempty(sprint(show, Forcing(forcing)))
end
