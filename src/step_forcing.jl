"""
Step-forcing container types and field-to-forcing conversion helpers.
"""

"""
    SnowpackStepForcing{NF}

Per-column forcing bundle consumed by [`step!`](@ref). It stores temperatures,
mass fluxes, optional prescribed surface-flux terms, and metadata needed by
the optional diurnal shortwave adjustment.
"""
struct SnowpackStepForcing{NF <: AbstractFloat}
    air_temperature::NF
    precipitation_rate::NF
    dt_days::NF
    snowfall_rate::NF
    rainfall_rate::NF
    shortwave_down::NF
    wind_speed::NF
    q_sw_net::NF
    q_lw_down::NF
    q_sh::NF
    q_lh::NF
    has_q_sw_net::Bool
    has_q_lw_down::Bool
    has_q_sh::Bool
    has_q_lh::Bool
    diurnal_shortwave::Bool
    latitude::NF
    day_of_year::NF
end

"""
    SnowpackStepFields

Batch forcing container that stores one forcing matrix per field for
multi-column, multi-time-step case execution.
"""
struct SnowpackStepFields{
        DT,
        AT,
        ST,
        RT,
        SWT,
        WST,
        QLWT,
        HQLWT,
        QSHT,
        HQSHT,
        QLHT,
        HQLHT,
    }
    dt_days::DT
    air_temperature::AT
    snowfall_rate::ST
    rainfall_rate::RT
    shortwave_down::SWT
    wind_speed::WST
    q_lw_down::QLWT
    has_q_lw_down::HQLWT
    q_sh::QSHT
    has_q_sh::HQSHT
    q_lh::QLHT
    has_q_lh::HQLHT
end

"""
    _assert_step_field_shape(name, field, field_shape)

Validate that `field` matches the `(ncol, ntime)` shape of the reference
forcing arrays. Throws an error with `name` when the shape is inconsistent.
"""
@inline function _assert_step_field_shape(name::AbstractString, field, field_shape)
    size(field) == field_shape || error("`$name` must match `air_temperature`.")
    return nothing
end

"""
    SnowpackStepFields(; ...)

Construct a batch forcing container from per-column, per-time-step arrays.
All matrix inputs must share the same shape as `air_temperature`, and
`dt_days` must either be scalar or provide one value per time step.
"""
function SnowpackStepFields(;
    dt_days,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    shortwave_down,
    wind_speed,
    q_lw_down,
    has_q_lw_down,
    q_sh,
    has_q_sh,
    q_lh,
    has_q_lh,
)
    field_shape = size(air_temperature)
    _assert_step_field_shape("snowfall_rate", snowfall_rate, field_shape)
    _assert_step_field_shape("rainfall_rate", rainfall_rate, field_shape)
    _assert_step_field_shape("shortwave_down", shortwave_down, field_shape)
    _assert_step_field_shape("wind_speed", wind_speed, field_shape)
    _assert_step_field_shape("q_lw_down", q_lw_down, field_shape)
    _assert_step_field_shape("has_q_lw_down", has_q_lw_down, field_shape)
    _assert_step_field_shape("q_sh", q_sh, field_shape)
    _assert_step_field_shape("has_q_sh", has_q_sh, field_shape)
    _assert_step_field_shape("q_lh", q_lh, field_shape)
    _assert_step_field_shape("has_q_lh", has_q_lh, field_shape)

    ntime = size(air_temperature, 2)
    if dt_days isa AbstractVector
        length(dt_days) == ntime || error("`dt_days` must have one entry per forcing time step.")
    else
        ntime == 1 || error("Scalar `dt_days` requires a single forcing time step.")
    end

    return SnowpackStepFields(
        dt_days,
        air_temperature,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        q_lw_down,
        has_q_lw_down,
        q_sh,
        has_q_sh,
        q_lh,
        has_q_lh,
    )
end

"""
    SnowpackStepFields(forcing)

Convert an equilibrium or script-level forcing container with matching field
names into a `SnowpackStepFields` instance.
"""
function SnowpackStepFields(forcing)
    return SnowpackStepFields(
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature,
        snowfall_rate=forcing.snowfall_rate,
        rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down,
        wind_speed=forcing.wind_speed,
        q_lw_down=forcing.q_lw_down,
        has_q_lw_down=forcing.has_q_lw_down,
        q_sh=forcing.q_sh,
        has_q_sh=forcing.has_q_sh,
        q_lh=forcing.q_lh,
        has_q_lh=forcing.has_q_lh,
    )
end

"""
    _step_time_count(fields)

Return the number of forcing time steps stored in `fields`.
"""
@inline _step_time_count(fields::SnowpackStepFields) = size(fields.air_temperature, 2)

"""
    _step_dt(dt_days, time_index)

Resolve the step duration in days for `time_index`, supporting both scalar and
vector-valued `dt_days` storage.
"""
@inline _step_dt(dt_days::Number, ::Int) = dt_days
@inline _step_dt(dt_days::AbstractVector, time_index::Int) = @inbounds dt_days[time_index]
@inline _step_dt(dt_days::CUDA.CuArray, time_index::Int) = CUDA.@allowscalar dt_days[time_index]

"""
    _step_forcing_from_fields(air_temperature, snowfall_rate, rainfall_rate, dt_days, shortwave_down, wind_speed, q_lw_down, has_q_lw_down, q_sh, has_q_sh, q_lh, has_q_lh)

Build a single-column `SnowpackStepForcing` from already-indexed forcing
values. The returned forcing disables optional fluxes that are not present and
sets precipitation rate to snowfall plus rainfall.
"""
@inline function _step_forcing_from_fields(
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    dt_days,
    shortwave_down,
    wind_speed,
    q_lw_down,
    has_q_lw_down::Bool,
    q_sh,
    has_q_sh::Bool,
    q_lh,
    has_q_lh::Bool,
)
    return SnowpackStepForcing(
        air_temperature,
        snowfall_rate + rainfall_rate,
        dt_days,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        zero(air_temperature),
        q_lw_down,
        q_sh,
        q_lh,
        false,
        has_q_lw_down,
        has_q_sh,
        has_q_lh,
        false,
        zero(air_temperature),
        zero(air_temperature),
    )
end

"""
    _step_forcing_from_fields(fields, idx, time_index)

Extract forcing values for column `idx` and time step `time_index` from
`fields` and package them as a `SnowpackStepForcing`.
"""
@inline function _step_forcing_from_fields(
    fields::SnowpackStepFields,
    idx::Int,
    time_index::Int,
)
    return _step_forcing_from_fields(
        @inbounds(fields.air_temperature[idx, time_index]),
        @inbounds(fields.snowfall_rate[idx, time_index]),
        @inbounds(fields.rainfall_rate[idx, time_index]),
        _step_dt(fields.dt_days, time_index),
        @inbounds(fields.shortwave_down[idx, time_index]),
        @inbounds(fields.wind_speed[idx, time_index]),
        @inbounds(fields.q_lw_down[idx, time_index]),
        @inbounds(fields.has_q_lw_down[idx, time_index]),
        @inbounds(fields.q_sh[idx, time_index]),
        @inbounds(fields.has_q_sh[idx, time_index]),
        @inbounds(fields.q_lh[idx, time_index]),
        @inbounds(fields.has_q_lh[idx, time_index]),
    )
end

"""
    SnowpackStepForcing(c, air_temperature, precipitation_rate, dt_days; ...)

Construct a single-step forcing bundle using physical constants `c` to choose
the model number type and default values. Optional turbulent and radiative
fluxes are tracked with explicit presence flags.
"""
function SnowpackStepForcing(
    c::SnowpackPhysicalConstants,
    air_temperature,
    precipitation_rate,
    dt_days;
    snowfall_rate=zero(air_temperature),
    rainfall_rate=zero(air_temperature),
    shortwave_down=oftype(air_temperature, 400.0),
    wind_speed=oftype(air_temperature, 10.0),
    q_sw_net=nothing,
    q_lw_down=nothing,
    q_sh=nothing,
    q_lh=nothing,
    diurnal_shortwave::Bool=false,
    latitude=zero(air_temperature),
    day_of_year=zero(air_temperature),
)
    NF = number_type(c)
    return SnowpackStepForcing(
        convert(NF, air_temperature),
        convert(NF, precipitation_rate),
        convert(NF, dt_days),
        convert(NF, snowfall_rate),
        convert(NF, rainfall_rate),
        convert(NF, shortwave_down),
        convert(NF, wind_speed),
        convert(NF, isnothing(q_sw_net) ? zero(NF) : q_sw_net),
        convert(NF, isnothing(q_lw_down) ? zero(NF) : q_lw_down),
        convert(NF, isnothing(q_sh) ? zero(NF) : q_sh),
        convert(NF, isnothing(q_lh) ? zero(NF) : q_lh),
        !isnothing(q_sw_net),
        !isnothing(q_lw_down),
        !isnothing(q_sh),
        !isnothing(q_lh),
        diurnal_shortwave,
        convert(NF, latitude),
        convert(NF, day_of_year),
    )
end

@adapt_structure SnowpackStepForcing
@adapt_structure SnowpackStepFields
