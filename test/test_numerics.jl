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
