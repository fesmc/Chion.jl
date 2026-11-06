#!/usr/bin/env python3
"""Model schematic of Chion for the model description paper.

(a) one model column: forcing, surface energy and mass fluxes, layer structure (near-surface
    geometric layers, snow/firn mass layers, thermal ice substrate) and internal processes;
(b) order of the processes in one (sub)daily step, as in `column_step_core!` (src/step.jl);
(c) implementation and coupling to an ice-sheet model.
Values are the BESSIModel defaults (src/models.jl, src/constants.jl).

Writes output/model_schematic/F0_model_schematic.{svg,pdf,png} (rsvg-convert).
    python3 examples/figures/model_schematic.py
"""
import os
import re
import subprocess

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "output", "model_schematic")
W, H = 1600, 1080
INK, MUTED = "#3a3a38", "#8a8a86"
BLUE, ORANGE, RED, GREEN, GREY = "#2a78d6", "#eda100", "#d6402a", "#1baf7a", "#8a8a86"
FONT = "DejaVu Sans"

svg = []
add = svg.append


def markup(s):
    """Escape XML and render T_s, T_air as subscripts."""
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return re.sub(r"([A-Za-z])_([a-z]+)", r'\1<tspan baseline-shift="sub" font-size="75%">\2</tspan>', s)


def text(x, y, s, size=21, anchor="start", weight="normal", fill=INK, rotate=None, style=""):
    tr = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    s = markup(s)
    add(f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" text-anchor="{anchor}" '
        f'font-weight="{weight}" fill="{fill}"{tr} {style}>{s}</text>')


def arrow(x1, y1, x2, y2, colour, width=3.5, dash=None, both=False):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = colour.lstrip("#")
    start = f' marker-start="url(#s{m})"' if both else ""
    add(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{colour}" stroke-width="{width}"{d} '
        f'marker-end="url(#e{m})"{start}/>')


def rect(x, y, w, h, fill, stroke="none", sw=1.5, rx=0, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')


def brace(x, y1, y2, label_lines, size=20):
    """Left-facing curly brace at x from y1 to y2 with the label to its left."""
    m = (y1 + y2) / 2
    add(f'<path d="M {x} {y1} q -10 0 -10 10 L {x-10} {m-8} q 0 8 -10 8 q 10 0 10 8 L {x-10} {y2-10} '
        f'q 0 10 10 10" fill="none" stroke="{INK}" stroke-width="1.6"/>')
    n = len(label_lines)
    for i, s in enumerate(label_lines):
        text(x - 30, m - (n - 1) * 12 + i * 24 + 7, s, size=size, anchor="end",
             weight="bold" if i == 0 else "normal")


add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
add("<defs>")
for c in (BLUE, ORANGE, RED, GREEN, GREY, INK):
    m = c.lstrip("#")
    add(f'<marker id="e{m}" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="4" markerHeight="4" '
        f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>')
    add(f'<marker id="s{m}" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="4" markerHeight="4" '
        f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>')
add('<linearGradient id="firn" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#eef2f7"/>'
    '<stop offset="1" stop-color="#b6c4d4"/></linearGradient>')
add('<linearGradient id="ice" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#b9d6ef"/>'
    '<stop offset="1" stop-color="#6fa3d4"/></linearGradient>')
add("</defs>")
rect(0, 0, W, H, "white")

# ============================================================== (a) column
text(20, 38, "(a)", size=24, weight="bold")
# forcing box
rect(20, 62, 330, 222, "#f6f6f4", stroke=MUTED, rx=8)
text(36, 92, "Forcing (daily, per column)", size=20, weight="bold")
for i, s in enumerate(["air temperature", "snowfall, rainfall", "shortwave radiation",
                       "relative humidity", "wind speed", "surface elevation, latitude"]):
    text(36, 124 + i * 28, "· " + s, size=20)

X0, X1, SURF = 400, 690, 440                       # column and surface
# surface fluxes: (x, direction, colour, dash, label)
fluxes = [(418, "down", BLUE, None, "Snowfall"), (452, "down", BLUE, "8 5", "Rainfall"),
          (492, "down", ORANGE, None, "Shortwave ↓"), (526, "up", ORANGE, "8 5", "Reflected (α)"),
          (566, "down", RED, None, "Longwave ↓"), (600, "up", RED, "8 5", "Longwave ↑"),
          (640, "both", GREEN, None, "Sensible heat"), (674, "both", GREEN, "8 5", "Latent heat, vapour")]
TAIL = 300
for x, d, c, dash, label in fluxes:
    if d == "down":
        arrow(x, TAIL, x, SURF - 6, c, dash=dash)
    elif d == "up":
        arrow(x, SURF - 4, x, TAIL + 6, c, dash=dash)
    else:
        arrow(x, TAIL + 6, x, SURF - 6, c, dash=dash, both=True)
    text(x + 7, TAIL - 10, label, size=19, rotate=-90)
# snowflakes and rain drops near the precipitation arrows
for y in (330, 380):
    text(398, y, "✻", size=18, fill=BLUE, anchor="middle")
# flux notes
rect(712, 120, 300, 212, "#f6f6f4", stroke=MUTED, rx=8)
for i, (s, c, b) in enumerate([("Albedo (dynamic)", ORANGE, True), ("dry 0.81, wet 0.70, ice 0.40;", INK, False),
                               ("ages in 20 d (cold), 2 d (melt)", INK, False),
                               ("Longwave ↓: cloud-proxy ε", RED, True),
                               ("Turbulent fluxes: SEMIX", GREEN, True), ("bulk scheme, stability-", INK, False),
                               ("dependent exchange", INK, False)]):
    text(724, 148 + i * 27, s, size=18, weight="bold" if b else "normal", fill=c if b else INK)

# layers
y = SURF
ns = [10, 14, 20, 30]
for h in ns:
    rect(X0, y, X1 - X0, h, "#fbfcfe", stroke=INK, sw=1.2)
    y += h
y_ns = y
for k in range(5):
    rect(X0, y, X1 - X0, 42, "url(#firn)", stroke=INK, sw=1.2)
    y += 42
y_firn_mid = y
text((X0 + X1) / 2, y + 22, "⋮", size=26, anchor="middle")
y += 30
rect(X0, y, X1 - X0, 36, "#aebdcf", stroke=INK, sw=1.2)
y_bottom_firn = y
y += 36
y_firn_end = y
y += 12
for h in [8, 12, 18, 28, 44]:
    rect(X0, y, X1 - X0, h, "url(#ice)", stroke=INK, sw=1.0)
    y += h
y_sub_end = y
rect(X0, y, X1 - X0, 46, "#6fa3d4", dash="6 6", stroke=INK, sw=1.0)
text((X0 + X1) / 2, y + 31, "ice sheet", size=20, anchor="middle", fill="white", weight="bold")
add(f'<line x1="{X0 - 40}" y1="{SURF}" x2="{X1 + 20}" y2="{SURF}" stroke="{INK}" stroke-width="2.5"/>')

# braces (left)
brace(X0 - 8, SURF, y_ns, ["Near-surface layers", "4 layers, ≤ 2, 5, 10, 30 cm"])
brace(X0 - 8, y_ns + 4, y_firn_end, ["Snow and firn", "mass layers (≤ 15 in total);", "split > 500 kg m⁻²,",
                                     "merge < 100 kg m⁻²"])
brace(X0 - 8, y_firn_end + 12, y_sub_end, ["Ice substrate", "5 thermal layers", "(top layer 5 cm)"])

# internal processes: percolation (blue), heat conduction (red)
for yy in (465, 505, 548, 590):
    add(f'<circle cx="470" cy="{yy}" r="5" fill="{BLUE}"/>')
arrow(470, 455, 470, y_bottom_firn - 4, BLUE, width=2.5, dash="5 5")
arrow(650, SURF + 8, 650, y_sub_end - 6, RED, width=2.5, dash="5 5", both=True)
arrow(X1, 640, X1 + 70, 640, BLUE, width=3.5)
text(X1 + 78, 647, "Runoff", size=20, weight="bold", fill=BLUE)
# bottom: full column -> ice
arrow(560, y_bottom_firn + 30, 560, y_firn_end + 22, INK, width=2.5)

# process labels (right)
labels = [(SURF + 30, "Surface temperature T_s;", RED), (SURF + 55, "implicit heat conduction", RED),
          (SURF + 95, "Melt at T_s = 0 °C", BLUE), (SURF + 135, "Percolation (bucket,", BLUE),
          (SURF + 160, "holding capacity)", BLUE), (SURF + 245, "Refreezing", BLUE),
          (SURF + 285, "Densification,", GREY), (SURF + 310, "wet compaction", GREY)]
for yy, s, c in labels:
    text(X1 + 20, yy, s, size=19, fill=INK)
text(X1 + 20, y_bottom_firn + 24, "Full column: bottom layer", size=19)
text(X1 + 20, y_bottom_firn + 48, "becomes ice (SMB)", size=19)
text(X1 + 20, y_sub_end - 30, "Bare ice: substrate", size=19)
text(X1 + 20, y_sub_end - 6, "melts (negative SMB)", size=19)

# subdaily note
rect(20, 985, 960, 72, "#fff6e0", stroke=ORANGE, rx=8)
text(36, 1013, "Subdaily resolution", size=20, weight="bold", fill=INK)
text(36, 1042, "≤ 8 shortwave substeps per day (SW↓ > 0, T_air > −8 °C); ±1 K diurnal T cycle", size=18)

# ============================================================== (b) daily step
BX = 1030
text(BX - 10, 38, "(b)", size=24, weight="bold")
text(BX + 42, 38, "One step of a column", size=21, weight="bold")
steps = [("Accumulation: snow, rain", BLUE), ("Near-surface remeshing", GREY), ("Albedo update", ORANGE),
         ("Densification", GREY), ("Energy balance: implicit T_s, T(z)", RED), ("Vapour exchange", GREEN),
         ("Melt", BLUE), ("Percolation and runoff", BLUE), ("Refreezing", BLUE),
         ("Remeshing, deep-layer splitting", GREY)]
for i, (s, c) in enumerate(steps):
    yy = 62 + i * 50
    rect(BX, yy, 450, 40, "white", stroke=c, sw=2.2, rx=6)
    rect(BX, yy, 40, 40, c, rx=6)
    text(BX + 20, yy + 28, str(i + 1), size=20, anchor="middle", fill="white", weight="bold")
    text(BX + 54, yy + 28, s, size=19)
yb = 62 + 10 * 50 - 10
add(f'<path d="M {BX + 470} 62 h 14 V {yb} h -14" fill="none" stroke="{ORANGE}" stroke-width="2.5"/>')
text(BX + 505, (62 + yb) / 2, "×1–8 substeps per day", size=19, anchor="middle", rotate=90, fill=INK)

# ============================================================== (c) implementation and coupling
CY = 640
text(BX - 10, CY, "(c)", size=24, weight="bold")
text(BX + 42, CY, "Implementation and coupling", size=21, weight="bold")
rect(BX, CY + 25, 140, 110, "#f6f6f4", stroke=MUTED, rx=8)
for i, s in enumerate(["Forcing", "MAR, ERA5,", "CMIP6 ΔT"]):
    text(BX + 70, CY + 58 + i * 26, s, size=19, anchor="middle", weight="bold" if i == 0 else "normal")
rect(BX + 175, CY + 25, 215, 110, "#eaf2fb", stroke=BLUE, sw=2, rx=8)
for i, s in enumerate(["Chion", "1 thread per column", "GPU or CPU,", "Float64 / Float32"]):
    text(BX + 282, CY + 52 + i * 25, s, size=18 if i else 19, anchor="middle", weight="bold" if i == 0 else "normal")
rect(BX + 425, CY + 25, 125, 110, "#f6f6f4", stroke=MUTED, rx=8)
for i, s in enumerate(["Yelmo", "ice-sheet", "dynamics"]):
    text(BX + 487, CY + 58 + i * 26, s, size=19, anchor="middle", weight="bold" if i == 0 else "normal")
arrow(BX + 142, CY + 80, BX + 172, CY + 80, INK, width=2.5)
arrow(BX + 392, CY + 62, BX + 422, CY + 62, INK, width=2.5)
add(f'<path d="M {BX + 487} {CY + 137} V {CY + 175} H {BX + 282} V {CY + 141}" fill="none" stroke="{INK}" '
    f'stroke-width="2.5" marker-end="url(#e{INK.lstrip("#")})"/>')
text(BX + 397, CY + 205, "annual: SMB, T_s → Yelmo;", size=18, anchor="middle")
text(BX + 397, CY + 230, "surface elevation → Chion", size=18, anchor="middle")
text(BX + 397, CY + 255, "(lapse-rate correction of T_air)", size=18, anchor="middle", fill=MUTED)
text(BX + 397, CY + 300, "Outputs: SMB, melt, runoff, refreezing,", size=18, anchor="middle")
text(BX + 397, CY + 325, "firn density and temperature, T_s", size=18, anchor="middle")

add("</svg>")

os.makedirs(OUT, exist_ok=True)
path = os.path.join(OUT, "F0_model_schematic.svg")
with open(path, "w") as f:
    f.write("\n".join(svg))
for fmt, extra in (("pdf", []), ("png", ["-z", "2"])):
    subprocess.run(["rsvg-convert", "-f", fmt, *extra, "-o", path[:-3] + fmt, path], check=True)
print("wrote", path[:-4] + ".{svg,pdf,png}")
