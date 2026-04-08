"""
Backend helpers for threaded and KernelAbstractions execution.
"""

"""
    kernelabstractions_available()

Report whether the package was built with KernelAbstractions support. This is
currently always `true` because the module unconditionally depends on
KernelAbstractions.
"""
@inline kernelabstractions_available() = true

"""
    cuda_available()

Return `true` when CUDA is functional in the current Julia session. This does
not allocate any device arrays; it only checks backend availability.
"""
@inline cuda_available() = CUDA.functional()

"""
    _ka_backend(array)

Return the KernelAbstractions backend associated with `array`. The result is
used to launch backend-specific kernels for CPU or GPU storage.
"""
@inline function _ka_backend(array)
    return KernelAbstractions.get_backend(array)
end

"""
    _wait_kernel(event)

Synchronize a KernelAbstractions launch event when one was returned. The
function is a no-op for `nothing` and always returns `nothing`.
"""
@inline function _wait_kernel(event)
    if !isnothing(event)
        KernelAbstractions.wait(event)
    end
    return nothing
end
