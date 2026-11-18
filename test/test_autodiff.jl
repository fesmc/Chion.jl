using Test
using Dates
using Enzyme
using Chion

# Melt-season forcing for a few columns, cold to warm, with snow on an ice substrate.
function _ad_forcing(ndays = 12)
    times = [DateTime(2001, 6, 20, 12) + Day(k - 1) for k in 1:ndays]
    return SnowpackForcing(dt_days = ones(ndays), ncol = 3, time_values = times,
        air_temperature_c = [-6.0 + 4.0c + 0.3k for c in 1:3, k in 1:ndays], snowfall_mm_day = 3.0,
        rainfall_mm_day = 0.5, shortwave_down = [220.0 + 15.0sin(k) for _ in 1:3, k in 1:ndays],
        wind_speed = 4.0, relative_humidity = 0.8, latitude_deg = 68.0, surface_height = 900.0)
end

# Snow at the start: the state after a spin-up of the same forcing.
function _ad_model(; kwargs...)
    model = BESSIModel(Grid(3; layering = Layering(Ntot = 8)); forcing = _ad_forcing(), kwargs...)
    integrate!(model, 1:12)
    return model
end

"""Melt plus runoff of column `idx` over the loaded year with the air temperature raised by `ΔT`,
stepping the state arrays `fields` with the kernels' `column_step!`."""
function melt_and_runoff!(fields, workspace, ΔT, model, idx)
    for time_index in 1:Chion.forcing_steps(model)
        step = Chion._step_forcing_at(model.forcing.fields, idx, time_index)
        step = Chion._with_step_forcing(step, (; air_temperature = step.air_temperature + ΔT))
        Chion.column_step!(fields, model.params, model.grid.layering, model.boundary, idx, step, workspace)
    end
    return fields.melt[idx] + fields.runoff[idx]
end

# Fresh copies of the model's state arrays and workspace, and zero tangents of the same shape.
_copies(model) = (map(copy, Chion.get_fields(model.state)), deepcopy(model.cache.workspace))
_zero_tangent(x::AbstractArray{<:AbstractFloat}) = zero(x)
_zero_tangent(x::AbstractArray) = copy(x)
_zero_tangent(x) = Chion._mapfields(_zero_tangent, x)

melt_and_runoff(ΔT, model, idx) = melt_and_runoff!(_copies(model)..., ΔT, model, idx)

function enzyme_slope(model, idx)
    fields, workspace = _copies(model)
    dfields = map(_zero_tangent, fields)
    dworkspace = _zero_tangent(workspace)
    return only(Enzyme.autodiff(Enzyme.Forward, melt_and_runoff!, Duplicated(fields, dfields),
        Duplicated(workspace, dworkspace), Duplicated(0.0, 1.0), Const(model), Const(idx)))
end

_central_difference(f, h) = (f(h) - f(-h)) / 2h

@testset "Enzyme forward mode through column_step!" begin
    for params in (Params(), Params(albedo_scheme = SemixAlbedo(), turbulent_flux_scheme = BESSITurbulentFluxes()))
        model = _ad_model(; params)
        for idx in 1:3
            slope = enzyme_slope(model, idx)
            fd = _central_difference(ΔT -> melt_and_runoff(ΔT, model, idx), 1e-6)
            @test isfinite(slope)
            # The objective is piecewise smooth (layer split/merge): AD gives the pathwise slope,
            # which finite differences resolve away from a switch.
            @test slope ≈ fd rtol = 1e-4 atol = 1e-6
        end
        @test melt_and_runoff(0.0, model, 3) > 0                       # the warm column melts
    end
end

using ForwardDiff

# A whole model at NF = Dual: the derivative of the melt and runoff of three columns over two
# passes of the forcing with respect to a parameter, through integrate! (CPU kernels).
function melt_and_runoff_total(value, name::Symbol)
    NF = typeof(value)
    params = Params(; NF, name => value)
    model = BESSIModel(Grid(3; NF, layering = Layering(Ntot = 8)); forcing = _ad_forcing(), params)
    integrate!(model, 1:12)
    integrate!(model, 1:12)
    return sum(model.state.melt) + sum(model.state.runoff)
end

@testset "ForwardDiff through a model at NF = Dual" begin
    for (name, value, h) in ((:alpha_dry, 0.81, 1e-7), (:alpha_ice, 0.40, 1e-7), (:refreezing_correction, 1.0, 1e-7))
        f(x) = melt_and_runoff_total(x, name)
        slope = ForwardDiff.derivative(f, value)
        fd = _central_difference(δ -> f(value + δ), h)
        @test isfinite(slope)
        @test slope ≈ fd rtol = 1e-4 atol = 1e-3
    end
    @test ForwardDiff.derivative(x -> melt_and_runoff_total(x, :alpha_dry), 0.81) < 0   # brighter snow melts less
end

@testset "Simulation with outputs at NF = Dual" begin
    NF = typeof(ForwardDiff.Dual(0.81, 1.0))
    model = BESSIModel(Grid(3; NF, layering = Layering(Ntot = 8)); forcing = _ad_forcing(),
        params = Params(; NF, alpha_dry = ForwardDiff.Dual(0.81, 1.0)))
    dir = mktempdir()
    sim = Simulation(model; verbose = false, output = OutputConfig(; dir, outputs = [
        Output("daily.nc"; variables = [:melt, :swe, :albedo, :density], frequency = :daily, kind = :snapshot),
        Output("profile.nc"; variables = :firn, frequency = :monthly, kind = :snapshot)]))
    run!(sim)
    @test finished(sim) && isfile(joinpath(dir, "daily.nc")) && isfile(joinpath(dir, "chion_run.toml"))
    @test ForwardDiff.partials(sum(model.state.melt))[1] != 0          # the state carries the derivative
    @test smbstats(model).mean_kg_m2.melt ≈ ForwardDiff.value(sum(model.state.melt)) / 3
end

@testset "SmoothKinks: a model runs and differentiates" begin
    f(x) = melt_and_runoff_total(x, :alpha_ice)
    smooth(x) = begin
        NF = typeof(x)
        params = Params(; NF, alpha_ice = x, kinks = SmoothKinks())
        model = BESSIModel(Grid(3; NF, layering = Layering(Ntot = 8)); forcing = _ad_forcing(), params)
        integrate!(model, 1:12); integrate!(model, 1:12)
        sum(model.state.melt) + sum(model.state.runoff)
    end
    slope = ForwardDiff.derivative(smooth, 0.40)
    @test isfinite(slope)
    @test slope ≈ _central_difference(δ -> smooth(0.40 + δ), 1e-7) rtol = 1e-4 atol = 1e-3
    @test smooth(0.40) ≈ f(0.40) rtol = 0.05                   # a small change of the result
end
