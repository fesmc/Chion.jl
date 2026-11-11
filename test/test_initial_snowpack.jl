using Test, Chion

const _CUMULATIVE_BESSI_FIELDS = (
    :mass_base, :smb_ice, :runoff, :melt, :refreezing, :vapor_mass, :sublimation,
    :latent_heat_flux_sum,
)

_bessi_array_fields(state) = (name => copy(getproperty(state, name)) for name in
    (:N, :mass, :mass_w, :density, :temperature, :Tsrf, :albedo, :snow_age_days,
     :thickness, :wet_mass, :bulk_density, :liquid_water, :ice_temperature,
     _CUMULATIVE_BESSI_FIELDS...))

_active(state, name, idx) = getproperty(state, name)[1:state.N[idx], idx]

function _remeshed(state, model, idx)
    copy_state = deepcopy(state)
    Chion._remesh_near_surface_layers!(copy_state.N, copy_state.mass, copy_state.mass_w,
        copy_state.density, copy_state.temperature, idx, model.Ntot,
        model.near_surface_layer_max_thicknesses_m, model.c)
    return copy_state
end

# Daily forcing for a cold and dry plateau column: -45 °C mean with a ±15 K
# seasonal cycle, 12 kg m⁻² yr⁻¹ of snowfall and 60% of the daily-mean
# top-of-atmosphere shortwave at 75° S.
function _cold_dry_forcing(; latitude_deg=-75.0, solar_constant=1361.0)
    days = 1:365
    air_temperature = [228.15 + 15cos(2π * (day - 15) / 365) for day in days]
    toa = [solar_constant / 2π * Chion._diurnal_shortwave_integral_terms(latitude_deg,
        Chion._solar_longitude_deg_from_calendar_day(Float64(day))).daylight_integral
        for day in days]
    forcing = SnowpackForcing(
        dt_days=ones(365),
        air_temperature=air_temperature,
        snowfall_rate=fill(12.0 / 365 / 86_400, 365),
        rainfall_rate=zeros(365),
        shortwave_down=0.6 .* toa,
        wind_speed=fill(5.0, 365),
        relative_humidity=fill(0.8, 365),
        latitude_deg=latitude_deg,
    )
    return forcing, sum(air_temperature) / 365
end

function _run_cold_dry_column(state_kwargs; years=100)
    forcing, mean_air_temperature = _cold_dry_forcing()
    model = BESSIModel(SnowpackGrid(1))
    sim = Simulation(model; forcing=forcing, state=initial_state(model; state_kwargs...), years=years)
    integrator = init_integrator(sim; io=devnull)
    annual = nothing
    for _ in 1:years
        annual = yearly_step!(integrator)
    end
    finalize!(integrator)
    update_diagnostics!(sim.ref)
    update_diagnostics!(sim.now)
    snowfall = years * 12.0
    Δstorage = sim.now.wet_mass[1] - sim.ref.wet_mass[1]
    return (
        model=model,
        storage_initial=sim.ref.wet_mass[1],
        storage=sim.now.wet_mass[1],
        accumulation_rate=Δstorage / years,
        mean_Tsrf_last_year=annual.mean_T_srf_K[1],
        mean_air_temperature=mean_air_temperature,
        albedo=sim.now.albedo[1],
        vapor_mass=sim.now.vapor_mass[1],
        budget_residual=snowfall + sim.now.vapor_mass[1] -
            (Δstorage + sim.now.runoff[1] + sim.now.smb_ice[1]),
    )
end

@testset "Initial snowpack (BESSI cold start)" begin
    @testset "no keywords: the snow-free state" begin
        model = BESSIModel(SnowpackGrid(3); Ntot=6)
        state = BESSIState(model)
        @test all(state.N .== 0)
        @test all(state.mass .== 0) && all(state.mass_w .== 0)
        @test all(state.density .== model.density_init)
        @test all(state.temperature .== model.temperature_init)
        @test all(state.ice_temperature .== model.temperature_init)
        @test all(state.Tsrf .== model.c.T0)
        @test all(state.albedo .== model.c.alpha_dry)
        for name in (:snow_age_days, :w_snow_max, :thickness, :wet_mass, :bulk_density,
                     :liquid_water, _CUMULATIVE_BESSI_FIELDS...)
            @test all(getproperty(state, name) .== 0)
        end
        for (keyword_state, label) in ((BESSIState(model; swe=nothing, T_init=nothing), :nothing),
                                       (initial_state(model), :initial_state))
            @test all(getproperty(keyword_state, name) == value
                      for (name, value) in _bessi_array_fields(state))
        end
    end

    @testset "swe 500 at 230 K: a layered firn column" begin
        model = BESSIModel(SnowpackGrid(2))
        state = BESSIState(model; swe=[500.0, 0.0], T_init=230.0)
        n = state.N[1]
        mass = _active(state, :mass, 1)
        thickness = mass ./ _active(state, :density, 1)
        caps = model.near_surface_layer_max_thicknesses_m

        @test sum(mass) ≈ 500.0 rtol=1e-14
        @test state.wet_mass[1] ≈ 500.0 rtol=1e-14
        @test 1 < n <= model.Ntot
        @test all(_active(state, :density, 1) .== model.density_init_firn)
        @test thickness[1:4] ≈ collect(caps)
        @test all(mass[5:n-1] .== model.mass_split)
        @test mass[n] >= model.mass_min
        @test sum(thickness) <= Chion.BESSI_REFERENCE_SNOW_DEPTH_M
        @test all(_active(state, :mass_w, 1) .== 0)
        # The accumulation-time remesh leaves the column as it is.
        remeshed = _remeshed(state, model, 1)
        @test remeshed.N[1] == n
        @test _active(remeshed, :mass, 1) ≈ mass rtol=1e-12

        @test all(_active(state, :temperature, 1) .== 230.0)
        @test state.Tsrf[1] == 230.0
        @test all(state.ice_temperature[:, 1] .== 230.0)
        @test state.albedo[1] == model.c.alpha_dry
        @test state.snow_age_days[1] == 0.0
        for name in _CUMULATIVE_BESSI_FIELDS
            @test all(getproperty(state, name) .== 0)
        end

        # swe = 0 with T_init: bare, substrate and surface at T_init.
        @test state.N[2] == 0
        @test all(state.mass[:, 2] .== 0)
        @test state.Tsrf[2] == 230.0
        @test all(state.ice_temperature[:, 2] .== 230.0)

        @test all(getproperty(initial_state(model; swe=[500.0, 0.0], T_init=230.0), name) == value
                  for (name, value) in _bessi_array_fields(state))
    end

    @testset "layering limits" begin
        model = BESSIModel(SnowpackGrid(4); Ntot=6)
        state = BESSIState(model; swe=[50.0, 200.0, 1850.0, 7000.0], T_init=[280.0, 250.0, 250.0, 250.0])
        # Thinner than the near-surface layers: a partial deepest layer.
        @test state.N[1] == 3
        @test sum(_active(state, :mass, 1)) ≈ 50.0
        @test all(_active(state, :temperature, 1) .== model.c.T0)
        @test state.Tsrf[1] == model.c.T0
        # A remainder below a fixed-thickness layer stays a layer.
        @test state.N[2] == 5
        @test state.mass[5, 2] ≈ 200.0 - sum(model.density_init_firn .* model.near_surface_layer_max_thicknesses_m)
        # A full column: the deepest layer takes the excess.
        @test state.N[4] == model.Ntot
        @test state.mass[model.Ntot, 4] ≈ 7000.0 - sum(state.mass[1:model.Ntot-1, 4])
        for idx in 1:4
            @test sum(_active(state, :mass, idx)) ≈ state.wet_mass[idx]
            remeshed = _remeshed(state, model, idx)
            @test remeshed.N[idx] == state.N[idx]
            @test _active(remeshed, :mass, idx) ≈ _active(state, :mass, idx) rtol=1e-12
        end

        # Without near-surface layers: mass_split layers, a small remainder merged up.
        coarse = BESSIModel(SnowpackGrid(2); Ntot=15, near_surface_layer_max_thicknesses_m=(Inf, Inf, Inf, Inf))
        coarse_state = BESSIState(coarse; swe=[1850.0, 2000.0], T_init=250.0)
        @test _active(coarse_state, :mass, 1) == [300.0, 300.0, 300.0, 300.0, 300.0, 350.0]
        @test _active(coarse_state, :mass, 2) == [300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 200.0]
    end

    @testset "invalid input" begin
        model = BESSIModel(SnowpackGrid(2))
        @test_throws ErrorException BESSIState(model; swe=100.0)
        @test_throws ErrorException BESSIState(model; swe=-1.0, T_init=250.0)
        @test_throws ErrorException BESSIState(model; swe=[100.0], T_init=250.0)
        @test_throws ErrorException BESSIState(model; swe=100.0, T_init=NaN)
        max_swe = model.density_init_firn * Chion.BESSI_REFERENCE_SNOW_DEPTH_M
        @test BESSIState(model; swe=max_swe, T_init=250.0).wet_mass[1] ≈ max_swe
        @test_throws ErrorException BESSIState(model; swe=1.01max_swe, T_init=250.0)
        @test_throws ErrorException BESSIModel(SnowpackGrid(1); density_init_firn=0.0)
        @test_throws ErrorException BESSIModel(SnowpackGrid(1); density_init_firn=1000.0)
    end

    @testset "cold-dry column accumulates from an initial firn column" begin
        firn = _run_cold_dry_column((swe=500.0, T_init=228.0))
        empty = _run_cold_dry_column(NamedTuple())
        @test firn.storage_initial == 500.0
        @test firn.accumulation_rate > 5.0
        @test abs(firn.mean_Tsrf_last_year - firn.mean_air_temperature) < 5.0
        @test firn.albedo ≈ firn.model.c.alpha_dry atol=0.02
        @test abs(firn.budget_residual) < 1e-8 * 1200
        # The same column started empty, for comparison (without the bare-surface
        # temperature fix it stays bare at Tsrf = T0).
        @test abs(empty.budget_residual) < 1e-8 * 1200
    end
end
