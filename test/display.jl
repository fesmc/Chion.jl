@testset "Display" begin
    forcing = SnowpackForcing(dt_days = [1.0, 1.0], ncol = 2, air_temperature_c = [-5.0, -4.0],
        snowfall_mm_day = 1.0, rainfall_mm_day = 0.0, shortwave_down = [100.0, 110.0], wind_speed = 3.0,
        latitude_deg = 70.0, time_values = [DateTime(2001, 1, 1, 12), DateTime(2001, 1, 2, 12)])
    grid = Grid(2)
    model = BESSIModel(grid; forcing, params = Params(albedo_scheme = AgingAlbedo(cold_timescale_days = 10)))
    sim = Simulation(model; verbose = false)
    plain(x) = sprint(show, MIME"text/plain"(), x)
    @test occursin("Simulation \"chion\" of a BESSIModel{Float64}", plain(sim))
    @test occursin("next: 2001 step 1", plain(sim))
    @test occursin("BESSIModel{Float64} on CPU", plain(model))
    @test occursin("year 2001 loaded", plain(model))
    @test occursin("Grid{Float64}: 2 columns", plain(grid))
    text = plain(model.params)
    @test occursin("Params{Float64}", text) && occursin("AgingAlbedo(cold_timescale_days = 10.0)", text)
    @test !occursin("alpha_dry", text)                      # numbers at their defaults are not listed
    @test occursin("alpha_dry", plain(Params(alpha_dry = 0.85)))
    @test sprint(show, SemixTurbulentFluxes()) == "SemixTurbulentFluxes()"
    @test sprint(show, DiurnalSubsteps(max_substeps = 4)) == "DiurnalSubsteps(max_substeps = 4)"
    @test occursin("PDDParams{Float32}", plain(PDDParams(NF = Float32)))
    @test !isempty(sprint(show, Forcing(forcing)))
end
