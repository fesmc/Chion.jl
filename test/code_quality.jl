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
