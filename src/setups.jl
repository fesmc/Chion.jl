#=
Named setups (as Laddie's reproductions): the simulations of published experiments, built in one
call and pinned by `test/reproductions.jl`.
=#

_setup_output(dir) = dir === nothing ? OutputConfig() : OutputConfig(; dir, outputs = [
    Output("monthly.nc"; variables = :monthly, frequency = :monthly),
    Output("budget_yearly.csv"; variables = :budget, frequency = :yearly, aggregate = :area_sum)])

"""
    SimulationGreenlandMARv1(path; years = 1940:2026, spinup = CycledSpinup(1940:1969; cycles = 5),
        surface_fluxes = :parameterized, files = Dict{Int, String}(), NF = Float64, backend = CPU(),
        params = ParamsGrISv1(; NF), dir = nothing, name = "greenland_mar_v1", verbose = true)

The Greenland run of the Chion paper: `BESSIModel` with the Greenland calibration on the ice columns
(MAR mask ≥ 50 %) of daily MAR v3.14.3 (10 km, ERA5), five spin-up cycles of 1940–1969, then
1940–2026. `path` is the pattern of the yearly MAR files with `{year}`; `files` replaces single
years (the paper used a repaired 2026 file). `surface_fluxes = :prescribed` (with
`params = ParamsGrISv1(; NF, albedo_scheme = PrescribedAlbedo())`) gives the prescribed-flux run.
With a run directory `dir`, writes `monthly.nc` and `budget_yearly.csv` there.
"""
function SimulationGreenlandMARv1(path::AbstractString; years = 1940:2026,
        spinup::AbstractSpinup = CycledSpinup(1940:1969; cycles = 5), surface_fluxes::Symbol = :parameterized,
        files = Dict{Int, String}(), NF::Type{<:Real} = Float64, backend = KernelAbstractions.CPU(),
        params = ParamsGrISv1(; NF), dir = nothing, name::AbstractString = "greenland_mar_v1", verbose::Bool = true)
    forcing = Forcing(path; format = :mar, surface_fluxes, files)
    model = BESSIModel(Grid(forcing; NF = number_type(params)); forcing, params, backend)
    return Simulation(model; years, spinup, output = _setup_output(dir), name, verbose)
end

"""
    SimulationAntarcticaMARv1(path, climatology; years = 1979:2026, spinup_years = 200,
        surface_fluxes = :parameterized, area_file = nothing, NF = Float64, backend = CPU(),
        params = ParamsGrISv1(; NF), dir = nothing, name = "antarctica_mar_v1", verbose = true)

The Antarctic run of the Chion paper: `BESSIModel` on the ice columns of daily MAR v3.14.3
(27.5 km, ERA5), spun up with `spinup_years` repetitions of the single-year file `climatology`
(the paper: 200 years of 1979), then 1979–2026. The paper used the Greenland calibration;
`params = ParamsAISv1(; NF)` gives the Antarctic one. `area_file` is a MAR file with the cell
areas. With a run directory `dir`, writes `monthly.nc` and `budget_yearly.csv` there.
"""
function SimulationAntarcticaMARv1(path::AbstractString, climatology::AbstractString; years = 1979:2026,
        spinup_years::Integer = 200, surface_fluxes::Symbol = :parameterized, area_file = nothing,
        NF::Type{<:Real} = Float64, backend = KernelAbstractions.CPU(), params = ParamsGrISv1(; NF),
        dir = nothing, name::AbstractString = "antarctica_mar_v1", verbose::Bool = true)
    forcing = Forcing(path; format = :mar, surface_fluxes, area_file)
    spinup_forcing = Forcing(climatology; format = :mar, surface_fluxes, area_file)
    model = BESSIModel(Grid(forcing; NF = number_type(params)); forcing, params, backend)
    spinup = spinup_years > 0 ? ClimatologySpinup(spinup_forcing; years = spinup_years) : NoSpinup()
    return Simulation(model; years, spinup, output = _setup_output(dir), name, verbose)
end
