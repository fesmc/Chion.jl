"""
Compute backends.

`RunOptions(backend=...)` selects a `ComputeBackend`; runtime preparation
dispatches on it. GPU support lives in the `ChionCUDAExt` package extension,
which is loaded by `using CUDA` and adds the `GPUBackend` methods of the hooks
below (`backend_functional`, `storage_type`, `_tune_ka_backend`,
`_multiprocessor_count`).
"""

abstract type ComputeBackend end
"""Multithreaded CPU execution (`backend=:threads`)."""
struct ThreadsBackend <: ComputeBackend end
"""CUDA GPU execution (`backend=:gpu`); requires `using CUDA`."""
struct GPUBackend <: ComputeBackend end

const _COMPUTE_BACKENDS = (threads=ThreadsBackend(), gpu=GPUBackend())

backend_name(::ThreadsBackend) = :threads
backend_name(::GPUBackend) = :gpu
Base.show(io::IO, backend::ComputeBackend) = print(io, backend_name(backend))

is_gpu(::ComputeBackend) = false
is_gpu(::GPUBackend) = true

"""Whether `backend` can run in the current session."""
backend_functional(::ThreadsBackend) = true
backend_functional(::ComputeBackend) = false

"""Array type holding model data on `backend`."""
storage_type(::ThreadsBackend) = Array
storage_type(backend::ComputeBackend) = _require_functional(backend)

_cuda_extension_loaded() = !isnothing(Base.get_extension(@__MODULE__, :ChionCUDAExt))

function _require_functional(backend::ComputeBackend)
    backend_functional(backend) && return backend
    is_gpu(backend) && !_cuda_extension_loaded() &&
        error("`backend=gpu` requires CUDA.jl: run `using CUDA` before setting up the simulation.")
    error("`backend=$(backend_name(backend))` requested, but it is not functional in the current environment.")
end

@inline cuda_available() = backend_functional(GPUBackend())
@inline gpu_storage_type() = storage_type(GPUBackend())

# KernelAbstractions backend for kernels launched on the storage of `array`.
@inline _ka_backend(array) = _tune_ka_backend(KernelAbstractions.get_backend(array))
@inline _tune_ka_backend(backend) = backend
@inline _tune_ka_backend(::KernelAbstractions.CPU) = KernelAbstractions.CPU(; static=true)

"""Multiprocessor count of the device behind `backend` (0 when unknown)."""
@inline _multiprocessor_count(::KernelAbstractions.Backend) = 0
