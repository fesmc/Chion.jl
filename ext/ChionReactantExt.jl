"""
Runs a `BESSIModel` through Reactant.jl: `to_backend(model, ReactantBackend())`. Loaded by
`using Reactant, CUDA`.

The state, workspace and forcing become Reactant arrays (the compact forcing matrices dense, the
layer fields transposed as on a GPU). The column step of one forcing step is compiled with the step
index as a traced number (the forcing is sliced at it with a dynamic slice, so the kernel runs at
local index 1), once per `_step!` call, and stepped from the host. A traced loop over the steps
(`@trace for`) does not compile yet (the kernel's custom call inside the while loop fails in
`raise_triton_custom_call`).
"""
module ChionReactantExt

using Chion
using Chion: BESSIModel, BESSIState, ColumnarStepWorkspace, KernelAbstractions
using Reactant
using Reactant: @compile, ConcreteRNumber, Ops
using CUDA

# Chion's `ReactantBackend()` is the backend of the transfer (Reactant's KernelAbstractions backend
# type is not available while this extension precompiles); models are recognised by their arrays.
Chion._reactant_ka_backend(::Nothing) = Chion.ReactantBackend()

Chion._to_device(::Chion.ReactantBackend, a::AbstractArray) = Reactant.to_rarray(Array(a))
Chion._forcing_to_device(::Chion.ReactantBackend, fields::NamedTuple) =
    map(x -> x isa AbstractArray ? Reactant.to_rarray(Array(x)) : x, fields)
Chion._forcing_to_device(::Chion.ReactantBackend, fields) =
    error("ReactantBackend supports BESSIModel only, not $(nameof(typeof(fields))).")
Chion._active_indices(::Union{Reactant.ConcretePJRTArray, Reactant.ConcreteIFRTArray}, indices) = Reactant.to_rarray(collect(Int, indices))

"""Kernels compiled with Reactant (models whose state arrays are Reactant arrays)."""
struct ReactantExecution end
Chion._execution(::Union{Reactant.ConcretePJRTArray, Reactant.ConcreteIFRTArray}) = ReactantExecution()

Chion._backend(::ReactantExecution, model) = Chion.ReactantBackend()

_at(x::Reactant.AnyTracedRArray{T, 1}, step) where {T} = Ops.dynamic_slice(x, [step], [1])
_at(x::Reactant.AnyTracedRArray{T, 2}, step) where {T} = Ops.dynamic_slice(x, [1, step], [size(x, 1), 1])
_at(x, step) = x

"""One forcing step `step` (traced) of the active columns: the forcing sliced at the step."""
function _step_at!(fields, params, layering, boundary, workspace, active_indices, forcing_fields, step)
    sliced = map(x -> _at(x, step), forcing_fields)
    Chion.launch!(Chion._step_columns_single_kernel!, fields.mass, length(active_indices),
        fields, params, layering, boundary, workspace, active_indices, sliced, 1)
    return nothing
end

function Chion._step!(::ReactantExecution, model::BESSIModel, time_start::Int, time_stop::Int)
    time_stop >= time_start || return nothing
    args = (Chion.get_fields(model.state), model.params, model.grid.layering, model.boundary,
        model.cache.workspace, model.cache.active_indices, model.forcing.fields)
    compiled = @compile _step_at!(args..., ConcreteRNumber(time_start))
    for time_index in time_start:time_stop
        compiled(args..., ConcreteRNumber(time_index))
    end
    return nothing
end

# Kernels outside the column step run on a host copy.
function Chion._update_diagnostics!(::ReactantExecution, model::BESSIModel)
    host = to_backend(model, Chion.CPU())
    Chion.update_diagnostics!(host)
    for name in (:thickness, :wet_mass, :bulk_density, :liquid_water)
        copyto!(getfield(model.cache, name), getfield(host.cache, name))
    end
    return model
end

Chion._reset_columns!(::ReactantExecution, model::BESSIModel, indices) =
    error("ReactantBackend cannot reset columns: use `set_active_mask!(…; reset_newly_inactive = false)`.")

end
