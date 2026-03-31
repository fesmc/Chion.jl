"""
Backend helpers for threaded and KernelAbstractions execution.
"""

@inline kernelabstractions_available() = true
@inline cuda_available() = CUDA.functional()

@inline function _ka_backend(array)
    return KernelAbstractions.get_backend(array)
end

@inline function _wait_kernel(event)
    if !isnothing(event)
        KernelAbstractions.wait(event)
    end
    return nothing
end

@inline function _fill_prefix!(buffer, value, n::Int)
    @inbounds for i in 1:n
        buffer[i] = value
    end
    return buffer
end

@inline function _copy_prefix!(dest, src, n::Int)
    @inbounds for i in 1:n
        dest[i] = src[i]
    end
    return dest
end
