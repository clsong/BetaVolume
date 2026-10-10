"""Generate the BetaVolume hex logo (man/figures/logo.svg and logo.png).

The sticker shows the idea of the method (Song et al. 2025): communities are
points in species space (axes A, B, C), the origin is always included, and
beta diversity is the volume of their convex hull. Here six of the seven
possible presence/absence compositions are present, so the hull is the unit
cube with one corner cut off; the dashed corner is the volume still missing to
reach the maximum beta diversity (d = 3). The red point is the origin.

Usage: python3 tools/make_logo.py   (requires numpy, scipy, cairosvg)
"""

import math
from pathlib import Path

import cairosvg
import numpy as np
from scipy.spatial import ConvexHull

OUT = Path(__file__).resolve().parent.parent / "man" / "figures"

# Hex geometry (pointy-top, standard 2 : sqrt(3) sticker proportions)
W, H = 518, 600
CX, CY = W / 2, H / 2
R = 292
hexagon = [(CX + R * math.cos(math.radians(a)), CY + R * math.sin(math.radians(a)))
           for a in range(-90, 270, 60)]

# Metacommunity: all compositions of 3 species except {A, C}, plus the origin
points = np.array([
    [0, 0, 0],  # origin
    [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [1, 1, 0], [0, 1, 1], [1, 1, 1],
], dtype=float)
missing = np.array([1, 0, 1], dtype=float)

# Orthographic camera looking from (+A, -B, +C): the origin and the cut-off
# facet {A, C, ABC} both face the viewer.
S = 158
ORIGIN = np.array([CX + 6, CY - 2])
view = np.array([1.0, -1.15, 0.85])
view /= np.linalg.norm(view)  # towards the viewer
right = np.cross([0.0, 0.0, 1.0], view)
right /= np.linalg.norm(right)
up = np.cross(view, right)


def proj(p):
    p = np.asarray(p, dtype=float) - 0.5  # centre the unit cube
    return ORIGIN + S * np.array([np.dot(p, right), -np.dot(p, up)])


hull = ConvexHull(points)

# Merge coplanar triangles into polygonal facets
facets = {}
for simplex, eq in zip(hull.simplices, hull.equations):
    key = tuple(np.round(eq, 6))
    facets.setdefault(key, set()).update(simplex.tolist())


def ordered(face_idx, normal):
    pts = points[list(face_idx)]
    c = pts.mean(axis=0)
    a = pts[0] - c
    a /= np.linalg.norm(a)
    b = np.cross(normal, a)
    ang = [math.atan2(np.dot(p - c, b), np.dot(p - c, a)) for p in pts]
    return [list(face_idx)[i] for i in np.argsort(ang)]


palette = {  # light, front-facing facets get lighter greens
    "top": "#9be7b8", "right": "#3fb98a", "left": "#1f8a70", "cut": "#c7f5dd",
}


def facet_color(n):
    if n[2] > 0.9:
        return palette["top"]
    if n[0] > 0.9:
        return palette["right"]
    if n[1] < -0.9:
        return palette["left"]
    return palette["cut"]


faces_svg = []
back_svg = []
for key, idx in facets.items():
    n = np.array(key[:3])
    poly = ordered(idx, n)
    d = " ".join(f"{proj(points[i])[0]:.1f},{proj(points[i])[1]:.1f}" for i in poly)
    if np.dot(n, view) > 1e-9:
        faces_svg.append(
            f'<polygon points="{d}" fill="{facet_color(n)}" fill-opacity="0.92" '
            f'stroke="#e9fff3" stroke-width="2.4" stroke-linejoin="round"/>')
    else:
        back_svg.append(
            f'<polygon points="{d}" fill="none" stroke="#e9fff3" stroke-opacity="0.25" '
            f'stroke-width="1.6" stroke-dasharray="5 6"/>')


def line(p, q, **kw):
    a, b = proj(p), proj(q)
    attrs = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in kw.items())
    return f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" {attrs}/>'


# Dashed outline of the missing corner (volume left to reach the maximum)
ghost = "".join(
    line(missing, q, stroke="#ffd36e", stroke_width=2.2, stroke_dasharray="6 6",
         stroke_linecap="round")
    for q in ([1, 0, 0], [0, 0, 1], [1, 1, 1]))

# Species axes from the origin
axes = []
for e, lab in zip(np.eye(3), ["A", "B", "C"]):
    tip = e * 1.28
    axes.append(line(np.zeros(3), tip, stroke="#e9fff3", stroke_opacity="0.55",
                     stroke_width=2, stroke_linecap="round"))
    t = proj(e * 1.42)
    axes.append(f'<text x="{t[0]:.1f}" y="{t[1] + 6:.1f}" font-family="DejaVu Sans" '
                f'font-size="20" font-style="italic" fill="#e9fff3" fill-opacity="0.75" '
                f'text-anchor="middle">{lab}</text>')

# Embedded communities (vertices) and the origin
dots = []
for p in points[1:]:
    q = proj(p)
    if np.allclose(p, [0, 1, 0]):  # behind the hull: draw faded
        dots.append(f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="6.5" fill="#5ec8ff" '
                    f'fill-opacity="0.35" stroke="#ffffff" stroke-opacity="0.35" stroke-width="1.5"/>')
        continue
    dots.append(f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="13" fill="#0b2e24" '
                f'fill-opacity="0.35"/>'
                f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="10" fill="#5ec8ff" '
                f'stroke="#ffffff" stroke-width="2.5"/>')
q = proj(points[0])
dots.append(f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="11" fill="#ff5a5f" '
            f'stroke="#ffffff" stroke-width="2.5"/>')

hex_pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in hexagon)
inner = [(CX + (R - 16) * math.cos(math.radians(a)), CY + (R - 16) * math.sin(math.radians(a)))
         for a in range(-90, 270, 60)]
inner_pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in inner)

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="#14523f"/>
      <stop offset="1" stop-color="#0a2a22"/>
    </linearGradient>
    <radialGradient id="glow" cx="0.5" cy="0.45" r="0.5">
      <stop offset="0" stop-color="#4fe0a0" stop-opacity="0.35"/>
      <stop offset="1" stop-color="#4fe0a0" stop-opacity="0"/>
    </radialGradient>
    <clipPath id="clip"><polygon points="{hex_pts}"/></clipPath>
  </defs>
  <polygon points="{hex_pts}" fill="url(#bg)"/>
  <g clip-path="url(#clip)">
    <circle cx="{CX}" cy="{CY - 20}" r="230" fill="url(#glow)"/>
    {"".join(back_svg)}
    {"".join(faces_svg)}
    {"".join(dots)}
    <text x="{CX}" y="{CY - 198}" font-family="DejaVu Serif" font-size="46" font-style="italic"
          fill="#ffd36e" text-anchor="middle">&#946;<tspan font-size="22" dy="9">vol</tspan></text>
    <text x="{CX}" y="{CY + 196}" font-family="DejaVu Sans" font-weight="bold" font-size="44"
          fill="#ffffff" text-anchor="middle" letter-spacing="0.5">BetaVolume</text>
  </g>
  <polygon points="{inner_pts}" fill="none" stroke="#9be7b8" stroke-opacity="0.35" stroke-width="2"/>
  <polygon points="{hex_pts}" fill="none" stroke="#9be7b8" stroke-width="14" stroke-linejoin="round"/>
</svg>
'''

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "logo.svg").write_text(svg)
cairosvg.svg2png(bytestring=svg.encode(), write_to=str(OUT / "logo.png"),
                 output_width=W, output_height=H)
print("wrote", OUT / "logo.svg", "and", OUT / "logo.png")
