"""Generate the ink-wash (水墨) variant of the BetaVolume hex logo.

Same geometry as tools/make_logo.py -- the convex hull of six presence/absence
compositions of three species plus the origin -- drawn in the manner of a
Chinese landscape painting: xuan-paper ground, ink-wash facets with mineral
blue-green (青绿) accents, distant mountains, a calligraphic inscription
"β多样性" (beta diversity) and a vermilion seal.

Usage: python3 tools/make_logo_ink.py   (requires numpy, scipy, pillow, cairosvg
and the "AR PL UKai" font, e.g. Debian/Ubuntu package fonts-arphic-ukai)
Writes man/figures/logo_ink.svg and man/figures/logo_ink.png.
"""

import base64
import io
import math
from pathlib import Path

import cairosvg
import numpy as np
from PIL import Image, ImageFilter
from scipy.spatial import ConvexHull

OUT = Path(__file__).resolve().parent.parent / "man" / "figures"
W, H = 518, 600
CX, CY = W / 2, H / 2
R = 292
rng = np.random.default_rng(7)

hexagon = [(CX + R * math.cos(math.radians(a)), CY + R * math.sin(math.radians(a)))
           for a in range(-90, 270, 60)]
hex_pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in hexagon)

# --- xuan paper texture (raster, embedded) ------------------------------------
grain = rng.normal(0, 1, (H, W))
grain = np.array(Image.fromarray(((grain + 4) * 30).clip(0, 255).astype(np.uint8))
                 .filter(ImageFilter.GaussianBlur(0.6)), dtype=float) / 255 - 0.5
fibres = np.zeros((H, W))
for _ in range(140):  # faint long fibres
    x0, y0 = rng.uniform(0, W), rng.uniform(0, H)
    ang = rng.uniform(0, math.pi)
    length = rng.uniform(15, 60)
    for t in np.linspace(0, length, int(length) * 2):
        x, y = int(x0 + t * math.cos(ang)), int(y0 + t * math.sin(ang))
        if 0 <= x < W and 0 <= y < H:
            fibres[y, x] += 1
fibres = np.array(Image.fromarray((fibres * 120).clip(0, 255).astype(np.uint8))
                  .filter(ImageFilter.GaussianBlur(0.8)), dtype=float) / 255
yy, xx = np.mgrid[0:H, 0:W]
vignette = ((xx - CX) ** 2 + (yy - CY) ** 2) / (R ** 2)
base = np.array([243, 234, 214], dtype=float)  # warm xuan paper
paper = base[None, None, :] + grain[..., None] * 10 - fibres[..., None] * 7 \
    - vignette[..., None] * np.array([18, 22, 30])
buf = io.BytesIO()
Image.fromarray(paper.clip(0, 255).astype(np.uint8)).save(buf, format="PNG")
paper_uri = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# --- distant mountains -----------------------------------------------------------
def ridge(y_base, amp, seed, n=9):
    r = np.random.default_rng(seed)
    xs = np.linspace(-20, W + 20, n)
    ys = y_base - amp * np.abs(r.normal(0.55, 0.35, n))
    ys[0] = ys[-1] = y_base - amp * 0.15
    pts = [(xs[0], H + 10)] + list(zip(xs, ys)) + [(xs[-1], H + 10)]
    d = f"M{pts[0][0]:.1f},{pts[0][1]:.1f} "
    # smooth with quadratic curves through midpoints
    for (x1, y1), (x2, y2) in zip(pts[1:-2], pts[2:-1]):
        d += f"Q{x1:.1f},{y1:.1f} {(x1 + x2) / 2:.1f},{(y1 + y2) / 2:.1f} "
    d += f"L{pts[-2][0]:.1f},{pts[-2][1]:.1f} L{pts[-1][0]:.1f},{pts[-1][1]:.1f} Z"
    return d


mountains = (
    f'<path d="{ridge(452, 95, 3)}" fill="url(#mist1)"/>'
    f'<path d="{ridge(492, 70, 11)}" fill="url(#mist2)"/>'
    f'<path d="{ridge(478, 26, 5)}" fill="url(#mist3)"/>'
)

# --- the geometric embedding --------------------------------------------------
points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                   [1, 1, 0], [0, 1, 1], [1, 1, 1]], dtype=float)
S = 128
ORIGIN = np.array([CX - 34, CY - 52])
view = np.array([1.0, -1.15, 0.85])
view /= np.linalg.norm(view)
right = np.cross([0.0, 0.0, 1.0], view)
right /= np.linalg.norm(right)
up = np.cross(view, right)


def proj(p):
    p = np.asarray(p, dtype=float) - 0.5
    return ORIGIN + S * np.array([np.dot(p, right), -np.dot(p, up)])


hull = ConvexHull(points)
facets = {}
for simplex, eq in zip(hull.simplices, hull.equations):
    facets.setdefault(tuple(np.round(eq, 6)), set()).update(simplex.tolist())


def ordered(idx, n):
    pts = points[list(idx)]
    c = pts.mean(axis=0)
    a = (pts[0] - c) / np.linalg.norm(pts[0] - c)
    b = np.cross(n, a)
    ang = [math.atan2(np.dot(p - c, b), np.dot(p - c, a)) for p in pts]
    return [list(idx)[i] for i in np.argsort(ang)]


def fill_for(n):
    if n[2] > 0.9:
        return "url(#washTop)"     # pale ink
    if n[0] > 0.9:
        return "url(#washRight)"   # azurite blue-green
    if n[1] < -0.9:
        return "url(#washLeft)"    # dense ink
    return "url(#washCut)"         # malachite green


faces, hidden = [], []
for key, idx in facets.items():
    n = np.array(key[:3])
    poly = ordered(idx, n)
    d = " ".join(f"{proj(points[i])[0]:.1f},{proj(points[i])[1]:.1f}" for i in poly)
    if np.dot(n, view) > 1e-9:
        faces.append(f'<polygon points="{d}" fill="{fill_for(n)}"/>')
        # brush outline: a dark stroke plus a dry-brush offset
        faces.append(f'<polygon points="{d}" fill="none" stroke="#1d1d1b" stroke-opacity="0.85" '
                     f'stroke-width="3.2" stroke-linejoin="round"/>')
        faces.append(f'<polygon points="{d}" fill="none" stroke="#1d1d1b" stroke-opacity="0.25" '
                     f'stroke-width="1.2" stroke-linejoin="round" transform="translate(1.6,1.2)"/>')
    else:
        hidden.append(f'<polygon points="{d}" fill="none" stroke="#1d1d1b" stroke-opacity="0.18" '
                      f'stroke-width="1.4" stroke-dasharray="4 6"/>')

dots = []
for p in points[1:]:
    q = proj(p)
    if np.allclose(p, [0, 1, 0]):
        dots.append(f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="5" fill="#1d1d1b" fill-opacity="0.25"/>')
        continue
    dots.append(f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="8.5" fill="#1d1d1b"/>'
                f'<circle cx="{q[0] - 2:.1f}" cy="{q[1] - 2:.1f}" r="2.4" fill="#f3ead8" fill-opacity="0.7"/>')
q = proj(points[0])  # origin in cinnabar
dots.append(f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="10" fill="#c8322b"/>'
            f'<circle cx="{q[0] - 2.5:.1f}" cy="{q[1] - 2.5:.1f}" r="2.6" fill="#f6d3c4" fill-opacity="0.8"/>')

# --- inscription and seal --------------------------------------------------------
IX, IY0, STEP = 398, 166, 44
inscription = (
    f'<text x="{IX}" y="{IY0}" font-family="DejaVu Serif" font-style="italic" font-size="38" '
    f'fill="#1d1d1b" text-anchor="middle">&#946;</text>'
)
for k, ch in enumerate("多样性"):
    inscription += (f'<text x="{IX}" y="{IY0 + STEP * (k + 1)}" font-family="AR PL UKai CN" '
                    f'font-size="40" fill="#1d1d1b" text-anchor="middle">{ch}</text>')

SX, SY, SS = IX - 22, IY0 + STEP * 3 + 22, 44
seal = (
    f'<g transform="rotate(-2 {SX + SS / 2} {SY + SS / 2})">'
    f'<rect x="{SX}" y="{SY}" width="{SS}" height="{SS}" rx="4" fill="#c8322b" fill-opacity="0.92"/>'
    f'<rect x="{SX + 3.5}" y="{SY + 3.5}" width="{SS - 7}" height="{SS - 7}" rx="2" fill="none" '
    f'stroke="#f3ead8" stroke-width="1.6" stroke-opacity="0.85"/>'
    f'<text x="{SX + SS / 2}" y="{SY + SS / 2 + 11}" font-family="AR PL UKai CN" font-size="31" '
    f'fill="#f3ead8" text-anchor="middle">宋</text></g>'
)

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     width="{W}" height="{H}" viewBox="0 0 {W} {H}">
  <defs>
    <clipPath id="clip"><polygon points="{hex_pts}"/></clipPath>
    <linearGradient id="mist1" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="#7d8e8a" stop-opacity="0.38"/>
      <stop offset="0.55" stop-color="#7d8e8a" stop-opacity="0.05"/>
    </linearGradient>
    <linearGradient id="mist2" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="#4f6b66" stop-opacity="0.45"/>
      <stop offset="0.6" stop-color="#4f6b66" stop-opacity="0.05"/>
    </linearGradient>
    <linearGradient id="mist3" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="#2a3633" stop-opacity="0.82"/>
      <stop offset="0.6" stop-color="#1d2523" stop-opacity="0.92"/>
    </linearGradient>
    <linearGradient id="washTop" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="#cfc6b1" stop-opacity="0.75"/>
      <stop offset="1" stop-color="#e9e0cb" stop-opacity="0.45"/>
    </linearGradient>
    <linearGradient id="washLeft" x1="0" y1="0" x2="0.4" y2="1">
      <stop offset="0" stop-color="#2b2b29" stop-opacity="0.92"/>
      <stop offset="1" stop-color="#57534b" stop-opacity="0.78"/>
    </linearGradient>
    <linearGradient id="washRight" x1="1" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="#2f6f86" stop-opacity="0.9"/>
      <stop offset="1" stop-color="#7aa6a8" stop-opacity="0.65"/>
    </linearGradient>
    <linearGradient id="washCut" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="#3f8f72" stop-opacity="0.85"/>
      <stop offset="1" stop-color="#a9cdb0" stop-opacity="0.55"/>
    </linearGradient>
  </defs>
  <g clip-path="url(#clip)">
    <image x="0" y="0" width="{W}" height="{H}" xlink:href="{paper_uri}"/>
    {mountains}
    {"".join(hidden)}
    {"".join(faces)}
    {"".join(dots)}
    {inscription}
    {seal}
    <text x="{CX}" y="{CY + 210}" font-family="DejaVu Serif" font-weight="bold" font-size="38"
          fill="#f3ead8" text-anchor="middle" letter-spacing="1">BetaVolume</text>
  </g>
  <polygon points="{hex_pts}" fill="none" stroke="#1d1d1b" stroke-width="12" stroke-linejoin="round"/>
  <polygon points="{" ".join(f"{CX + (R - 13) * math.cos(math.radians(a)):.1f},{CY + (R - 13) * math.sin(math.radians(a)):.1f}" for a in range(-90, 270, 60))}"
           fill="none" stroke="#c8322b" stroke-width="2" stroke-opacity="0.8"/>
</svg>
'''

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "logo_ink.svg").write_text(svg)
cairosvg.svg2png(bytestring=svg.encode(), write_to=str(OUT / "logo_ink.png"),
                 output_width=W, output_height=H)
print("wrote", OUT / "logo_ink.svg", "and", OUT / "logo_ink.png")
