#=
Compute backends. Models live on a KernelAbstractions backend: `CPU()` (threads), or
`CUDABackend()` from CUDA.jl, whose support is the `ChionCUDAExt` package extension (`using CUDA`).
The names `:threads` and `:gpu` (keywords, TOML) map to them through `_ka_backend_from`.
=#

"""Name of a compute backend in keywords and TOML files: `:threads` or `:gpu`."""
abstract type ComputeBackend end
"""Multithreaded CPU execution (`:threads`)."""
struct ThreadsBackend <: ComputeBackend end
"""CUDA GPU execution (`:gpu`); requires `using CUDA`."""
struct GPUBackend <: ComputeBackend end

_compute_backend(name::Symbol) = name === :threads ? ThreadsBackend() : name === :gpu ? GPUBackend() :
    throw(ArgumentError("Unknown backend `:$name`; use `:threads` or `:gpu`."))

"""Whether `backend` can run in the current session (the CUDA extension adds the GPU)."""
backend_functional(::ThreadsBackend) = true
backend_functional(::ComputeBackend) = false

_cuda_extension_loaded() = !isnothing(Base.get_extension(@__MODULE__, :ChionCUDAExt))

"""
    _ka_backend_from(backend)

The KernelAbstractions backend of `backend`: a KernelAbstractions backend as it is, or the name
`:threads` / `:gpu`.
"""
_ka_backend_from(backend::KernelAbstractions.Backend) = backend
_ka_backend_from(name::Union{Symbol, AbstractString}) = _ka_backend_from(_compute_backend(Symbol(name)))
_ka_backend_from(::ThreadsBackend) = KernelAbstractions.CPU()
_ka_backend_from(::GPUBackend) = _cuda_ka_backend(nothing)

# The CUDA extension adds the method for `::Nothing`.
_cuda_ka_backend(::Any) = error("`backend = :gpu` requires CUDA.jl: run `using CUDA` first.")

"""
    cuda_available() -> Bool

Whether the GPU backend can be used: CUDA.jl is loaded (`using CUDA`, which activates the
`ChionCUDAExt` extension) and a functional GPU is present.
"""
@inline cuda_available() = backend_functional(GPUBackend())

# KernelAbstractions backend for kernels launched on the storage of `array`.
@inline _ka_backend(array) = _tune_ka_backend(KernelAbstractions.get_backend(array))
@inline _tune_ka_backend(backend) = backend
@inline _tune_ka_backend(::KernelAbstractions.CPU) = KernelAbstractions.CPU(; static=true)

"""Multiprocessor count of the device behind `backend` (0 when unknown)."""
@inline _multiprocessor_count(::KernelAbstractions.Backend) = 0

"""Workgroup size chosen by KernelAbstractions (the `launch!` default)."""
struct DefaultWorkgroup end

@inline _instantiate(kernel!, backend, ::DefaultWorkgroup) = kernel!(backend)
@inline _instantiate(kernel!, backend, workgroupsize) = kernel!(backend, workgroupsize)

"""
    launch!(kernel!, reference, ndrange, args...; workgroupsize = DefaultWorkgroup())

Launch the KernelAbstractions kernel `kernel!` with `args` over `ndrange` work items on the
backend of the array `reference` (the CPU backend in its static-threads form, see
`_tune_ka_backend`). An empty range launches nothing. The one place where kernels are
instantiated, as Laddie's `launch!`.
"""
@inline function launch!(kernel!, reference, ndrange, args...; workgroupsize = DefaultWorkgroup())
    prod(ndrange) == 0 && return nothing
    backend = _ka_backend(reference)
    _instantiate(kernel!, backend, workgroupsize)(args...; ndrange)
    return nothing
end
