"""
Liquid-water percolation following the original BESSI-style routine.
"""

"""
    go_percolation!(
        snowman::Vector{Float64},
        lwmass::Vector{Float64},
        rho_snow::Vector{Float64},
        rho_i::Float64,
        rho_w::Float64;
        max_lwc::Float64 = 0.05,
        rho_i_tol::Float64 = 10.0,
    ) -> Float64

Translate the Fortran `go_percolation` logic to Julia.

- `snowman`: layer snow mass [kg m^-2]
- `lwmass`: layer liquid-water mass [kg m^-2]
- `rho_snow`: layer snow density [kg m^-3]
- returns `runoff` produced by percolation in this call [kg m^-2]
"""
function go_percolation!(
    snowman::AbstractVector{Float64},
    lwmass::AbstractVector{Float64},
    rho_snow::AbstractVector{Float64},
    rho_i::Float64,
    rho_w::Float64;
    max_lwc::Float64 = 0.05,
    rho_i_tol::Float64 = 10.0,
)
    n_snowlayer = length(snowman)
    @assert length(lwmass) == n_snowlayer
    @assert length(rho_snow) == n_snowlayer

    runoff = 0.0
    ii = 1
    while ii <= n_snowlayer
        if snowman[ii] > 0.0
            if rho_snow[ii] > rho_i - rho_i_tol
                # Very dense snow: push all liquid water downward.
                percolating = lwmass[ii]
                lwmass[ii] = 0.0
                if ii < n_snowlayer
                    if snowman[ii + 1] > 0.0
                        lwmass[ii + 1] += percolating
                    else
                        runoff += percolating
                    end
                else
                    runoff += percolating
                end
            else
                lwc = lwmass[ii] / snowman[ii] / rho_w / (1.0 / rho_snow[ii] - 1.0 / rho_i)
                if lwc > max_lwc
                    percolating = (lwc - max_lwc) * rho_w * snowman[ii] * (1.0 / rho_snow[ii] - 1.0 / rho_i)
                    lwmass[ii] -= percolating
                    if ii < n_snowlayer
                        if snowman[ii + 1] > 0.0
                            lwmass[ii + 1] += percolating
                        else
                            runoff += percolating
                        end
                    else
                        runoff += percolating
                    end
                end
            end
        else
            # Keep Fortran control flow: abort loop at first empty box.
            ii = n_snowlayer
        end
        ii += 1
    end

    return runoff
end

"""
    go_percolation!(column::SnowpackColumn; max_lwc=0.05, rho_i_tol=10.0) -> Float64

Apply percolation to active layers of a `SnowpackColumn`.
Returned runoff is also added to `column.runoff`.
"""
function go_percolation!(
    column::SnowpackColumn;
    max_lwc::Float64 = 0.05,
    rho_i_tol::Float64 = 10.0,
)
    n = column.N
    if n <= 0
        return 0.0
    end

    @views runoff = go_percolation!(
        column.mass[1:n],
        column.mass_w[1:n],
        column.density[1:n],
        column.c.rho_i,
        column.c.rho_w;
        max_lwc=max_lwc,
        rho_i_tol=rho_i_tol,
    )
    column.runoff += runoff
    return runoff
end
