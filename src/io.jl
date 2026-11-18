#=
Metadata of the state variables of the bulk models and of the NetCDF time axis, used by the
output recorders (`outputs.jl`).
=#

struct OutputMetadata
    name::String
    long_name::String
    units::String
end

const NETCDF_METADATA = Dict{Symbol, OutputMetadata}(
    :thickness => OutputMetadata("thickness", "Snow thickness", "m"),
    :wet_mass => OutputMetadata("wet_mass", "Snow wet mass", "mmWE"),
    :bulk_density => OutputMetadata("bulk_density", "Bulk snow density", "kg m-3"),
    :liquid_water => OutputMetadata("liquid_water", "Liquid water mass", "kg m-2"),
    :mass_base => OutputMetadata("mass_base", "Firn mass exported to the ice model", "mmWE"),
    :smb_ice => OutputMetadata("smb_ice", "Net mass forcing to the ice sheet", "mmWE"),
    :surface_smb => OutputMetadata("surface_smb", "Surface mass balance", "mmWE"),
    :runoff => OutputMetadata("runoff", "Cumulative runoff", "mmWE"),
    :melt => OutputMetadata("melt", "Cumulative melt", "mmWE"),
    :refreezing => OutputMetadata("refreezing", "Cumulative refreezing", "mmWE"),
    :sublimation => OutputMetadata("sublimation", "Cumulative sublimation", "mmWE"),
    :latent_heat_flux => OutputMetadata("latent_heat_flux", "Monthly mean turbulent latent heat flux", "W m-2"),
    :latent_heat_flux_sum => OutputMetadata("latent_heat_flux_sum", "Integrated turbulent latent heat flux", "W m-2"),
    :shortwave_down => OutputMetadata("shortwave_down", "Monthly mean shortwave down", "W m-2"),
    :shortwave_net => OutputMetadata("shortwave_net", "Monthly mean absorbed shortwave", "W m-2"),
    :longwave_down => OutputMetadata("longwave_down", "Monthly mean longwave down", "W m-2"),
    :longwave_net => OutputMetadata("longwave_net", "Monthly mean net longwave", "W m-2"),
    :sensible_heat_flux => OutputMetadata("sensible_heat_flux", "Monthly mean sensible heat flux", "W m-2"),
    :rain_heat_flux => OutputMetadata("rain_heat_flux", "Monthly mean heat of rain", "W m-2"),
    :shortwave_down_sum => OutputMetadata("shortwave_down_sum", "Integrated shortwave down", "W m-2 day"),
    :shortwave_net_sum => OutputMetadata("shortwave_net_sum", "Integrated absorbed shortwave", "W m-2 day"),
    :longwave_down_sum => OutputMetadata("longwave_down_sum", "Integrated longwave down", "W m-2 day"),
    :longwave_net_sum => OutputMetadata("longwave_net_sum", "Integrated net longwave", "W m-2 day"),
    :sensible_heat_flux_sum => OutputMetadata("sensible_heat_flux_sum", "Integrated sensible heat flux", "W m-2 day"),
    :rain_heat_flux_sum => OutputMetadata("rain_heat_flux_sum", "Integrated heat of rain", "W m-2 day"),
    :Tsrf => OutputMetadata("Tsrf", "Surface temperature", "K"),
    :albedo => OutputMetadata("albedo", "Surface albedo", "1"),
    :snow_age_days => OutputMetadata("snow_age_days", "Time since the latest snowfall event", "day"),
    :N => OutputMetadata("N", "Number of active snow layers", "1"),
    :mass => OutputMetadata("mass", "Layer snow mass", "kg m-2"),
    :mass_w => OutputMetadata("mass_w", "Layer liquid-water mass", "kg m-2"),
    :density => OutputMetadata("density", "Layer density", "kg m-3"),
    :temperature => OutputMetadata("temperature", "Layer temperature", "K"),
    :snowpack_swe => OutputMetadata("snowpack_swe", "Snowpack water equivalent", "mmWE"),
    :pdd_sum => OutputMetadata("pdd_sum", "Cumulative positive degree days", "degC day"),
    :H_snow => OutputMetadata("H_snow", "ITM snowpack water equivalent", "mmWE"),
    :smb => OutputMetadata("smb", "ITM total surface mass balance rate", "mmWE day-1"),
    :smbi => OutputMetadata("smbi", "ITM ice-facing mass balance rate", "mmWE day-1"),
    :melt_net => OutputMetadata("melt_net", "ITM net melt rate", "mmWE day-1"),
    :smb_cum => OutputMetadata("smb_cum", "ITM cumulative surface mass balance", "mmWE"),
    :melt_cum => OutputMetadata("melt_cum", "ITM cumulative melt", "mmWE"),
    :runoff_cum => OutputMetadata("runoff_cum", "ITM cumulative runoff", "mmWE"),
    :refreezing_cum => OutputMetadata("refreezing_cum", "ITM cumulative refreezing", "mmWE"),
    :alb_s => OutputMetadata("alb_s", "ITM surface albedo", "1"),
)

# Monthly mean surface energy fluxes (W m⁻²) and the state sums (W m⁻² day)
# they are taken from.
const MONTHLY_FLUX_SUMS = (
    :latent_heat_flux => :latent_heat_flux_sum,
    :shortwave_down => :shortwave_down_sum,
    :shortwave_net => :shortwave_net_sum,
    :longwave_down => :longwave_down_sum,
    :longwave_net => :longwave_net_sum,
    :sensible_heat_flux => :sensible_heat_flux_sum,
    :rain_heat_flux => :rain_heat_flux_sum,
)
const MONTHLY_OUTPUT_VARS = (:surface_smb, :smb_ice, :runoff, :melt, :refreezing, :sublimation,
    first.(MONTHLY_FLUX_SUMS)..., :albedo)
const PDD_OUTPUT_VARS = (:snowpack_swe, :smb_ice, :runoff, :pdd_sum)
const ITM_OUTPUT_VARS = (:H_snow, :alb_s, :smb, :smbi, :melt, :runoff, :refreezing, :Tsrf, :melt_net, :smb_cum, :smb_ice, :melt_cum, :runoff_cum, :refreezing_cum)

"""State variables of a bulk model that outputs can record."""
output_variables(::PDDModel) = PDD_OUTPUT_VARS
output_variables(::ITMModel) = ITM_OUTPUT_VARS

@inline _host_vector(data::Vector{Float64}; copy_array::Bool=false) = copy_array ? copy(data) : data
@inline _host_vector(data; copy_array::Bool=false) = Float64.(Array(data))

const NETCDF_TIME_EPOCH = DateTime(1970, 1, 1)
const NETCDF_TIME_UNITS = "days since 1970-01-01 00:00:00"

_netcdf_time_days(value::DateTime) = Dates.value(value - NETCDF_TIME_EPOCH) / 86_400_000
