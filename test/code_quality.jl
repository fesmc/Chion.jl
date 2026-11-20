using Aqua

@testset "Code quality (Aqua)" begin
    Aqua.test_all(Chion)
end

_source_files() = [joinpath(root, f) for (root, _, fs) in walkdir(joinpath(pkgdir(Chion), "src")) for f in fs
    if endswith(f, ".jl")]
_code(f) = join(filter(l -> !startswith(lstrip(l), "#"), readlines(f)), "\n")
_relative(f) = relpath(f, joinpath(pkgdir(Chion), "src"))

# Launches go through `launch!`.
@testset "Code conventions" begin
    for f in _source_files()
        endswith(f, "backends.jl") && continue
        @testset "$(_relative(f))" begin
            # Kernels are instantiated only in `launch!` (src/backends.jl).
            @test !occursin(r"_kernel!\(\s*(_ka_backend|backend)\b", _code(f))
        end
    end
end

# Choices are types selected by dispatch: no scheme flags or symbol comparisons in the model code,
# and the removed containers stay removed.
@testset "Schemes by dispatch" begin
    for f in _source_files()
        text = read(f, String)
        @testset "$(_relative(f))" begin
            @test !occursin(r"::UInt8\b", text)
            @test !occursin(r"_uses_semix|SEB_SEMIX|SEB_BESSI|LONGWAVE_GRAYBODY|SEMIX_ALBEDO_", text)
            @test !occursin(r"isa\s+(Constant|Dynamic|Prescribed|Aging|Semix)Albedo\b", text)
            @test !occursin(r"SnowpackPhysicalConstants|BESSIParameters|SnowpackGrid|RunOptions|SimulationIntegrator", text)
        end
    end
end

# Every parameter and scheme field documents its unit and default, as the parameter pages show them.
@testset "Parameter documentation" begin
    for P in (Params, PDDParams, ITMParams, Layering, ConstantFreshSnowDensity, ParameterizedFreshSnowDensity,
            AgingAlbedo, SemixAlbedo, SemixTurbulentFluxes, BESSITurbulentFluxes, CloudProxyLongwave,
            DiurnalSubsteps, DiurnalTemperatureCycle)
        docs = string(Base.Docs.doc(P))
        for name in fieldnames(P)
            @test occursin(string(name), docs)
        end
    end
end
