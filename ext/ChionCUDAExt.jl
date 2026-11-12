"""CUDA GPU backend for Chion, loaded by `using CUDA`."""
module ChionCUDAExt

using Chion
using CUDA

Chion.backend_functional(::Chion.GPUBackend) = CUDA.functional()
Chion.storage_type(backend::Chion.GPUBackend) =
    CUDA.functional() ? CUDA.CuArray : Chion._require_functional(backend)

# Inline the whole column step into the CUDA kernel: without it, every non-inlined call
# passes large structs (e.g. the 0.5 kB physical constants) through local memory, which
# spilled 2.4 kB per thread; inlining cuts the spills to 0.95 kB and is ~1.3x faster.
@inline Chion._tune_ka_backend(::CUDA.CUDABackend) = CUDA.CUDABackend(; always_inline=true)

const MULTIPROCESSORS = Ref(0)
@inline function Chion._multiprocessor_count(::CUDA.CUDABackend)
    MULTIPROCESSORS[] > 0 && return MULTIPROCESSORS[]
    MULTIPROCESSORS[] = CUDA.attribute(CUDA.device(), CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)
    return MULTIPROCESSORS[]
end

end
