using Aqua

@testset "Code quality (Aqua)" begin
    Aqua.test_all(Chion)
end

# Choices are types selected by dispatch; launches go through `launch!`.
@testset "Code conventions" begin
    src = joinpath(pkgdir(Chion), "src")
    files = [joinpath(root, f) for (root, _, fs) in walkdir(src) for f in fs if endswith(f, ".jl")]
    code(f) = join(filter(l -> !startswith(lstrip(l), "#"), readlines(f)), "\n")
    # Kernels are instantiated only in `launch!` (src/backends.jl).
    for f in files
        endswith(f, "backends.jl") && continue
        @testset "$(relpath(f, src))" begin
            @test !occursin(r"_kernel!\(\s*(_ka_backend|backend)\b", code(f))
        end
    end
end

# No scheme is selected by a flag or a symbol comparison inside the model code.
@testset "Schemes by dispatch" begin
    src = joinpath(pkgdir(Chion), "src")
    for f in (joinpath(root, f) for (root, _, fs) in walkdir(src) for f in fs if endswith(f, ".jl"))
        text = read(f, String)
        @testset "$(relpath(f, src))" begin
            @test !occursin(r"::UInt8\b", text)
            @test !occursin(r"_uses_semix|SEB_SEMIX|SEB_BESSI|LONGWAVE_GRAYBODY|SEMIX_ALBEDO_", text)
            @test !occursin(r"isa\s+(Constant|Dynamic|Prescribed|Aging|Semix)Albedo\b", text)
        end
    end
end
