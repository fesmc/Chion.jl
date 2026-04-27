module ChionExampleScriptHelpers

using Chion

export arg_value, has_flag
export _albedo_scheme, _densification_scheme, _fresh_snow_scheme

function arg_value(args::Vector{String}, name::String, default::String="")
    prefix = "--" * name * "="
    for arg in args
        startswith(arg, prefix) && return arg[length(prefix)+1:end]
    end
    return default
end

has_flag(args::Vector{String}, name::String) = any(==("--" * name), args)

_albedo_scheme(name) =
    lowercase(strip(String(name))) == "constant" ? Chion.ConstantAlbedo() : Chion.DynamicAlbedo()

_densification_scheme(name) =
    lowercase(strip(String(name))) == "htessel" ? Chion.HTESSELDensification() : Chion.BESSIDensification()

_fresh_snow_scheme(name) =
    lowercase(strip(String(name))) == "parameterized" ?
        Chion.ParameterizedFreshSnowDensity() :
        Chion.ConstantFreshSnowDensity()

end
