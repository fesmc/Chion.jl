#=
Backend-generic batched tridiagonal solver for layer-by-column arrays.

The coefficient convention matches the energy solver:
  lower[k, col] multiplies row k + 1 and layer k,
  diagonal[k, col] multiplies row k and layer k,
  upper[k, col] multiplies row k and layer k + 1.
=#

struct BatchedTridiagonalSolver{LT,DT,UT,ST,NT,AT,AM}
    lower::LT
    diagonal::DT
    upper::UT
    scratch::ST
    active_layers::NT
    active_indices::AT
    active_mask::AM
end

function BatchedTridiagonalSolver(;
    lower_diagonal,
    diagonal,
    upper_diagonal,
    scratch,
    active_layers,
    active_indices=axes(diagonal, 2),
    active_mask=nothing,
)
    return BatchedTridiagonalSolver(
        lower_diagonal,
        diagonal,
        upper_diagonal,
        scratch,
        active_layers,
        active_indices,
        active_mask,
    )
end

@inline _tridiagonal_workgroupsize(backend) =
    backend isa KernelAbstractions.CPU ? 512 : 256

@inline _tridiagonal_column_is_active(::Nothing, ::Int) = true
@inline _tridiagonal_column_is_active(active_mask, idx::Int) = @inbounds active_mask[idx]

@inline function _solve_tridiagonal_column!(
    solution,
    lower_diagonal,
    main_diagonal,
    upper_diagonal,
    right_hand_side,
    scratch,
    idx::Int,
    n::Int,
)
    n <= 0 && return solution

    @inbounds begin
        β = _get_layer(main_diagonal, 1, idx)
        _set_layer!(solution, 1, idx, _get_layer(right_hand_side, 1, idx) / β)

        for layer_index in 2:n
            previous_upper = _get_layer(upper_diagonal, layer_index - 1, idx)
            lower = _get_layer(lower_diagonal, layer_index - 1, idx)
            _set_layer!(scratch, layer_index, idx, previous_upper / β)

            β = _get_layer(main_diagonal, layer_index, idx) -
                lower * _get_layer(scratch, layer_index, idx)

            _set_layer!(
                solution,
                layer_index,
                idx,
                (
                    _get_layer(right_hand_side, layer_index, idx) -
                    lower * _get_layer(solution, layer_index - 1, idx)
                ) / β,
            )
        end

        for layer_index in (n - 1):-1:1
            _set_layer!(
                solution,
                layer_index,
                idx,
                _get_layer(solution, layer_index, idx) -
                _get_layer(scratch, layer_index + 1, idx) *
                _get_layer(solution, layer_index + 1, idx),
            )
        end
    end

    return solution
end

@kernel function _solve_batched_tridiagonal_kernel!(
    solution,
    lower_diagonal,
    main_diagonal,
    upper_diagonal,
    right_hand_side,
    scratch,
    active_layers,
    active_indices,
    active_mask,
    max_layers::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        idx = active_indices[active_idx]
        if _tridiagonal_column_is_active(active_mask, idx)
            n = min(_n_active(active_layers, idx), max_layers)
            _solve_tridiagonal_column!(
                solution,
                lower_diagonal,
                main_diagonal,
                upper_diagonal,
                right_hand_side,
                scratch,
                idx,
                n,
            )
        end
    end
end

function solve_batched_tridiagonal!(
    solution,
    solver::BatchedTridiagonalSolver,
    right_hand_side;
    max_layers::Int=size(solution, 1),
)
    backend = _ka_backend(solution)
    kernel! = _solve_batched_tridiagonal_kernel!(backend, _tridiagonal_workgroupsize(backend))
    _wait_kernel(kernel!(
        solution,
        solver.lower,
        solver.diagonal,
        solver.upper,
        right_hand_side,
        solver.scratch,
        solver.active_layers,
        solver.active_indices,
        solver.active_mask,
        max_layers,
        ndrange=length(solver.active_indices),
    ))
    return solution
end

function _solve_batched_tridiagonal_prefix!(
    solution,
    lower_diagonal,
    main_diagonal,
    upper_diagonal,
    right_hand_side,
    scratch,
    active_layers,
    active_indices;
    max_layers::Int=size(solution, 1),
    active_mask=nothing,
)
    solver = BatchedTridiagonalSolver(;
        lower_diagonal,
        diagonal=main_diagonal,
        upper_diagonal,
        scratch,
        active_layers,
        active_indices,
        active_mask,
    )
    return solve_batched_tridiagonal!(solution, solver, right_hand_side; max_layers)
end
