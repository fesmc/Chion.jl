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
end
