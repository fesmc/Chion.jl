# Shared figure style and helpers for the paper figures (same conventions as
# evaluation_figures.jl): DejaVu Sans, no titles (panel letter above the axes, the panel's
# subject inside), fixed categorical colours, Crameri maps for fields.
using Plots, Printf, Statistics, Dates

const C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
const INK, MUTED = "#3a3a38", "#8a8a86"
const SSP_COLOUR = Dict("ssp126" => C1, "ssp245" => C4, "ssp585" => C2)
const SSP_LABEL = Dict("ssp126" => "SSP1-2.6", "ssp245" => "SSP2-4.5", "ssp585" => "SSP5-8.5")
const MODEL_COLOUR = Dict("chion" => C1, "pdd" => C2, "itm" => C3, "mar" => INK)
const MODEL_LABEL = Dict("chion" => "Chion", "pdd" => "PDD", "itm" => "ITM", "mar" => "MAR")
Plots.default(fontfamily="DejaVu Sans", titlefontsize=9, guidefontsize=8, tickfontsize=7,
    legendfontsize=7, framestyle=:box, grid=false, linewidth=1.5, foreground_color_legend=nothing,
    background_color_legend=nothing, titlelocation=:left, dpi=300)

panel_tag(i) = "($(Char('a' + i - 1)))"
"""Text inside the axes at a corner (row 2, 3, … stack below). Call after the limits are final."""
function note!(p, label; loc=:topleft, size=8, row=1)
    (x0, x1), (y0, y1) = Plots.xlims(p), Plots.ylims(p)
    left, top = loc in (:topleft, :bottomleft), loc in (:topleft, :topright)
    annotate!(p, x0 + (left ? 0.03 : 0.97) * (x1 - x0), y0 + (top ? 1.05 - 0.08 * row : 0.08 * row - 0.05) * (y1 - y0),
        text(label, size, INK, left ? :left : :right, top ? :top : :bottom))
    return p
end
"""Panel letter above the axes and the panel's subject inside them."""
tag!(p, i, label=""; loc=:topleft) = (title!(p, panel_tag(i)); isempty(label) ? p : note!(p, label; loc))
"""A legend-only strip for multi-panel figures; entries are (label, colour, keywords)."""
function legend_strip(entries; columns=length(entries))
    p = plot(; framestyle=:none, legend=:top, legend_column=columns, xlims=(0, 1), ylims=(0, 1))
    for (label, colour, kw) in entries; plot!(p, [NaN], [NaN]; c=colour, label, kw...); end
    return p
end
with_headroom!(p, f=0.14) = (yl = Plots.ylims(p); ylims!(p, (yl[1], yl[2] + f * (yl[2] - yl[1]))); p)
function save_figure(fig, dir, name)
    mkpath(dir)
    savefig(fig, joinpath(dir, "$name.pdf")); savefig(fig, joinpath(dir, "$name.png"))
end

"""Statistics collected while plotting and written to stats.csv."""
const STATS = NamedTuple[]
stat!(figure, quantity, value, unit) = push!(STATS, (; figure, quantity, value, unit))
function write_stats(dir)
    open(joinpath(dir, "stats.csv"), "w") do io
        println(io, "figure,quantity,value,unit")
        for s in STATS; println(io, join((s.figure, '"' * s.quantity * '"', s.value, '"' * s.unit * '"'), ',')); end
    end
end

"""CSV as Dict(column => Vector{Float64}); non-numeric cells become NaN (text kept in :_text)."""
function read_table(path)
    lines = filter(l -> !isempty(l) && !startswith(l, '#'), readlines(path))
    header = strip.(split(lines[1], ','))
    cells = [split(l, ',') for l in lines[2:end]]
    t = Dict{Symbol,Any}(Symbol(h) => [(v = tryparse(Float64, c[i]); isnothing(v) ? NaN : v) for c in cells]
                         for (i, h) in enumerate(header))
    t[:_text] = Dict(Symbol(h) => [String(c[i]) for c in cells] for (i, h) in enumerate(header))
    return t
end
trend_per_decade(x, y) = 10 * sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2)
running(v, w=11) = [(s = filter(isfinite, v[max(1, i - w ÷ 2):min(end, i + w ÷ 2)]); isempty(s) ? NaN : mean(s)) for i in eachindex(v)]
nice(x) = (e = 10.0^floor(log10(x)); m = x / e; (m <= 1 ? 1 : m <= 2 ? 2 : m <= 2.5 ? 2.5 : m <= 5 ? 5 : 10) * e)
