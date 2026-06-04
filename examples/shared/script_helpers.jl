module ChionExampleScriptHelpers

using Chion

export arg_value, env_value, has_flag
export _albedo_scheme, _densification_scheme, _fresh_snow_scheme

function arg_value(args::Vector{String}, name::String, default::String="")
    prefix = "--" * name * "="
    for arg in args
        startswith(arg, prefix) && return arg[length(prefix)+1:end]
    end
    return default
end

function env_value(name::String, default::String="")
    value = get(ENV, name, "")
    return isempty(value) ? default : value
end

has_flag(args::Vector{String}, name::String) = any(==("--" * name), args)

function _albedo_scheme(name)
    normalized = lowercase(strip(String(name)))
    normalized == "constant" && return Chion.ConstantAlbedo()
    normalized == "dynamic" && return Chion.DynamicAlbedo()
    normalized == "prescribed" && return Chion.PrescribedAlbedo()
    error("Unsupported albedo scheme '$(name)'. Use constant, dynamic, or prescribed.")
end

_densification_scheme(name) =
    lowercase(strip(String(name))) == "htessel" ? Chion.HTESSELDensification() : Chion.BESSIDensification()

_fresh_snow_scheme(name) =
    lowercase(strip(String(name))) == "parameterized" ?
        Chion.ParameterizedFreshSnowDensity() :
        Chion.ConstantFreshSnowDensity()

end
