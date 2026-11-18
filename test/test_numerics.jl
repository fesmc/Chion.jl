using Test
using KernelAbstractions: @kernel, @index
using Chion

@kernel function _test_fill_kernel!(x, value)
    i = @index(Global)
    @inbounds x[i] = value
end

@testset "Numerical guards" begin
    for T in (Float64, Float32)
        # _safe_pow is the plain power wherever that is defined (bitwise), and 0 at 0.
        for x in T[0, 1e-30, 1e-6, 0.3, 0.5, 0.999, 1], p in T[0.125, 1 / 2.75, 1, 3]
            @test Chion._safe_pow(x, p) === (x == 0 ? zero(T) : x^p)
        end
        @test Chion._safe_pow(-zero(T), T(0.5)) === zero(T)
        @test Chion._safe_positive(zero(T)) === T(Chion.EPS_TINY)
        @test Chion._safe_positive(T(NaN)) === T(Chion.EPS_TINY)
        @test Chion._safe_positive(T(2)) === T(2)
        @test Chion._safe_nonnegative(T(-1)) === zero(T)
        @test Chion._safe_nonnegative(T(NaN)) === zero(T)
        @test Chion._safe_nonnegative(T(3)) === T(3)
    end
end

@testset "launch!" begin
    x = zeros(5)
    Chion.launch!(_test_fill_kernel!, x, length(x), x, 2.0)
    @test all(==(2.0), x)
    Chion.launch!(_test_fill_kernel!, x, 3, x, 4.0; workgroupsize = 2)
    @test x == [4.0, 4.0, 4.0, 2.0, 2.0]
    @test Chion.launch!(_test_fill_kernel!, x, 0, x, 9.0) === nothing   # empty range: nothing launched
    @test x == [4.0, 4.0, 4.0, 2.0, 2.0]
end

@testset "Params" begin
    p = Params()
    @test isbitstype(typeof(p))
    @test !any(T -> T <: Integer, fieldtypes(typeof(p)))
    @test p.albedo_scheme === DynamicAlbedo() && p.diurnal_cycle == DiurnalSubsteps()
    # Every number and every float of a scheme in the number type.
    p32 = Params(NF = Float32, turbulent_flux_scheme = BESSITurbulentFluxes(D_sh = 8),
        albedo_scheme = SemixAlbedo(snow_bands = WarrenWiscombeSnowAlbedo(), frac_vu = 0.5))
    @test p32.alpha_dry === 0.81f0
    @test p32.turbulent_flux_scheme === BESSITurbulentFluxes{Float32}(8f0, 1f0)
    @test p32.albedo_scheme.frac_vu === 0.5f0 && p32.albedo_scheme.snow_bands === WarrenWiscombeSnowAlbedo()
    @test p32.diurnal_cycle.max_substeps === 8 && p32.diurnal_cycle.min_air_temperature === Float32(273.15 - 8.0)
    @test p32.diurnal_cycle.temperature_cycle.amplitude === 1f0
    # Integer literals become floats; unknown fields and bad values throw.
    @test AgingAlbedo(cold_timescale_days = 15) === AgingAlbedo(15.0, 2.0)
    @test DiurnalSubsteps(max_substeps = 4, threshold = 1) isa DiurnalSubsteps{Float64}
    @test Layering(mass_max = 600).mass_max === 600.0
    # All-integer coefficients become the number type of the parameters (integer fields stay Int).
    p_int = Params(NF = Float32, albedo_scheme = AgingAlbedo(cold_timescale_days = 15, melting_timescale_days = 2),
        diurnal_cycle = DiurnalSubsteps(max_substeps = 4, threshold = 0, min_air_temperature = 265,
            temperature_cycle = DiurnalTemperatureCycle(amplitude = 1, amplitude_gradient = 0, reference_height = 0,
                amplitude_max = 2)))
    @test p_int.albedo_scheme === AgingAlbedo{Float32}(15f0, 2f0)
    @test p_int.diurnal_cycle.max_substeps === 4 && p_int.diurnal_cycle.threshold === 0f0
    @test p_int.diurnal_cycle.temperature_cycle.amplitude_max === 2f0
    @test Grid(1; layering = Layering(mass_max = 600, mass_split = 300, mass_min = 100,
        near_surface_layer_max_thicknesses_m = (1, 1, 1, 1), ice_substrate_top_thickness_m = 1)).mass_max === 600.0
    @test_throws ArgumentError Params(rho_s = 300.0)
    @test_throws ArgumentError Params(diurnal_cycle = DiurnalSubsteps(max_substeps = 30))
    @test_throws ArgumentError Params(turbulent_flux_scheme = SemixTurbulentFluxes(karman = 0))
    @test_throws ArgumentError Params(longwave_scheme = CloudProxyLongwave(night_cloud_fraction = 2))
    @test_throws ArgumentError Params(albedo_scheme = AgingAlbedo(cold_timescale_days = 0))
    @test_throws ArgumentError Grid(2; layering = Layering(mass_split = 600.0))
    # Copy with changes; named setups.
    @test Params(p; alpha_dry = 0.85).alpha_dry == 0.85
    @test ParamsGrISv1() == Params()
    ais = ParamsAISv1()
    @test ais.alpha_dry == 0.87 && ais.turbulent_flux_scheme.stable_coefficient == 20.0
    @test ais.longwave_scheme.night_cloud_fraction == 0.9 && ais.longwave_scheme.clear_sky_emissivity == 0.7314
    @test PDDParams(method = PISMPDD()).method === PISMPDD()
    @test_throws ArgumentError PDDParams(ddf_snow = -1)
    @test_throws ArgumentError ITMParams(alb_ice = 2)
end

@testset "Grid" begin
    grid = Grid(3; NF = Float32, layering = Layering(Ntot = 20, ice_substrate_layers = 0))
    @test grid.Ntot == 20 && grid.ice_substrate_layers == 0 && grid.mass_max === 500f0
    @test grid.x === nothing && !Chion.has_spatial_coords(grid)
    @test Grid(grid; layering = Layering()).Ntot == 15
    spatial = Grid(2; x = [0.0, 1.0], y = [0.0], js = [1, 1], is = [1, 2])
    @test spatial.mask == ones(1, 2) && Chion.has_spatial_coords(spatial)
    @test_throws ArgumentError Grid(2; x = [0.0, 1.0])
    @test_throws ArgumentError Grid(0)
end
