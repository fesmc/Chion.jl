# The Reactant backend against the CPU. Own environment (Reactant is heavy; its tests must not run
# with --check-bounds=yes, as Pkg.test does):
#   julia --project=test/reactant -e 'using Pkg; Pkg.instantiate()'
#   julia --project=test/reactant test/reactant/runtests.jl
using Test, Dates
using Reactant, CUDA
using Chion

function _forcing(ndays = 10)
    times = [DateTime(2001, 6, 20, 12) + Day(k - 1) for k in 1:ndays]
    return SnowpackForcing(dt_days = ones(ndays), ncol = 5, time_values = times,
        air_temperature_c = [-6.0 + 2.5c + 0.3k for c in 1:5, k in 1:ndays], snowfall_mm_day = 3.0,
        rainfall_mm_day = 0.5, shortwave_down = [220.0 + 15.0sin(k) for _ in 1:5, k in 1:ndays],
        wind_speed = 4.0, relative_humidity = 0.8, latitude_deg = 68.0, surface_height = 900.0)
end

_model(; kwargs...) = BESSIModel(Grid(5; layering = Layering(Ntot = 8)); forcing = _forcing(), kwargs...)

@testset "Reactant backend" begin
    cpu = _model()
    integrate!(cpu, 1:10)
    reactant = to_backend(_model(), ReactantBackend())
    @test reactant.state.mass.parent isa Reactant.AnyConcreteRArray
    integrate!(reactant, 1:10)
    host = cpu_state(reactant)
    for name in (:mass, :density, :temperature, :melt, :runoff, :refreezing, :Tsrf, :albedo, :ice_temperature)
        @test getfield(host, name) ≈ getfield(cpu.state, name) rtol = 1e-10 atol = 1e-9
    end
    @test host.N == cpu.state.N

    # A Simulation (no outputs) and the diagnostics on the Reactant model.
    sim = Simulation(_model(; backend = ReactantBackend()); verbose = false)
    run!(sim)
    @test finished(sim)
    @test Array(sim.model.thickness) ≈ (update_diagnostics!(cpu); cpu.thickness) rtol = 1e-10
    back = to_backend(sim.model, CPU())
    @test back.state.mass ≈ cpu.state.mass rtol = 1e-10 atol = 1e-9
end
