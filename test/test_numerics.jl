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

@testset "Scheme selection by dispatch" begin
    S = Chion
    @test S._scheme(S.AbstractAlbedoScheme, :semix) === S.SemixAlbedo()
    @test S._scheme(S.AbstractAlbedoScheme, "aging") === S.AgingAlbedo()
    @test S._scheme(S.AbstractAlbedoScheme, S.ConstantAlbedo()) === S.ConstantAlbedo()
    @test S._scheme(S.AbstractLongwaveScheme, :graybody) === S.GraybodyLongwave()
    @test S._scheme(S.AbstractSurfaceEnergyBalanceScheme, :bessi) === S.BESSISurfaceEnergyBalance()
    @test S._scheme(S.AbstractTurbulentFluxScheme, :semix) === S.SemixTurbulentFluxes()
    @test S._scheme(S.AbstractSemixSnowAlbedo, :warren_wiscombe) === S.WarrenWiscombeSnowAlbedo()
    @test S._scheme(S.AbstractPDDMethod, :pism) === S.PISMPDD()
    err = try S._scheme(S.AbstractAlbedoScheme, :nonsense); catch e; e; end
    @test err isa ArgumentError
    @test occursin("Unsupported albedo scheme 'nonsense'", err.msg)
    @test occursin(":constant, :dynamic, :prescribed, :aging, :semix", err.msg)
    @test_throws ArgumentError S._scheme(S.AbstractAlbedoScheme, S.GraybodyLongwave())   # wrong kind
    @test_throws ArgumentError S._scheme(S.AbstractPDDMethod, :calov_greve)
    # names round-trip through the table
    for (kind, _, schemes) in S._SCHEME_TABLE, (name, scheme) in pairs(schemes)
        @test S._scheme(kind, name) === scheme
        @test S._scheme_name(scheme) === name
    end
    # every choice is a type parameter: no scheme flags left in the constants
    c = S.SnowpackPhysicalConstants(; seb_scheme = :bessi, longwave_scheme = :graybody,
        turbulent_flux_scheme = S.BESSITurbulentFluxes(), semix_snow_albedo = :warren_wiscombe)
    @test !any(T -> T <: Integer, fieldtypes(typeof(c)))
    @test S.surface_energy_balance_scheme(c) === S.BESSISurfaceEnergyBalance()
    @test S.longwave_scheme(c) === S.GraybodyLongwave()
    @test S.turbulent_flux_scheme(c) === S.BESSITurbulentFluxes()
    @test S.semix_snow_albedo(c) === S.WarrenWiscombeSnowAlbedo()
    @test isbitstype(typeof(c))
    @test_throws ArgumentError S.SnowpackPhysicalConstants(; albedo_scheme = :aging, alpha_wet = 0.9, alpha_dry = 0.8)
end
