#!/usr/bin/env python3
"""
3D extruded US state map — actual state polygon shapes, height = mean_elevation * group_mult.
107W east/west split, W/E = 0.3321.  Gold diamond = mass centroid near Laramie, WY.
"""
import json, math
from urllib.request import urlopen
import numpy as np
import plotly.graph_objects as go
from scipy.spatial import Delaunay

# ── State data ──────────────────────────────────────────────────────────────
STATES = [
    ("Alabama",        32.80,  -86.80,  135767,  150),
    ("Alaska",         64.20, -153.40, 1723337,  580),
    ("Arizona",        34.29, -111.66,  295234, 1247),
    ("Arkansas",       34.89,  -92.44,  137732,  200),
    ("California",     37.17, -119.47,  423968,  884),
    ("Colorado",       38.99, -105.55,  269601, 2073),
    ("Connecticut",    41.62,  -72.73,   14357,   91),
    ("Delaware",       39.00,  -75.50,    6446,   18),
    ("Florida",        28.63,  -81.52,  170312,   30),
    ("Georgia",        32.64,  -83.44,  153910,  183),
    ("Hawaii",         20.25, -156.33,   28313,  916),
    ("Idaho",          44.38, -114.61,  216443, 1528),
    ("Illinois",       40.04,  -89.20,  149995,  179),
    ("Indiana",        40.27,  -86.13,   94327,  228),
    ("Iowa",           42.08,  -93.50,  145746,  335),
    ("Kansas",         38.53,  -98.38,  213100,  610),
    ("Kentucky",       37.52,  -85.30,  104656,  230),
    ("Louisiana",      31.07,  -91.96,  135659,   30),
    ("Maine",          45.37,  -68.97,   91633,  182),
    ("Maryland",       39.06,  -76.80,   32131,   91),
    ("Massachusetts",  42.26,  -71.81,   27336,  122),
    ("Michigan",       44.35,  -85.41,   96714,  274),
    ("Minnesota",      46.31,  -93.09,  225163,  366),
    ("Mississippi",    32.74,  -89.67,  125438,   91),
    ("Missouri",       38.46,  -92.29,  180540,  229),
    ("Montana",        47.03, -110.45,  380831, 1012),
    ("Nebraska",       41.49,  -99.68,  200330,  792),
    ("Nevada",         39.33, -116.63,  286382, 1676),
    ("New Hampshire",  43.68,  -71.58,   24214,  305),
    ("New Jersey",     40.14,  -74.67,   22591,   58),
    ("New Mexico",     34.84, -106.25,  314917, 1682),
    ("New York",       42.94,  -75.52,  141297,  305),
    ("North Carolina", 35.54,  -79.38,  139391,  244),
    ("North Dakota",   47.44, -100.47,  183123,  595),
    ("Ohio",           40.42,  -82.79,  116099,  262),
    ("Oklahoma",       35.59,  -96.93,  181036,  396),
    ("Oregon",         43.93, -120.56,  254806, 1097),
    ("Pennsylvania",   40.59,  -77.21,  119280,  284),
    ("Rhode Island",   41.68,  -71.52,    4001,   61),
    ("South Carolina", 33.92,  -81.10,   82933,  104),
    ("South Dakota",   44.37, -100.35,  199730,  731),
    ("Tennessee",      35.86,  -86.35,  109153,  281),
    ("Texas",          31.49,  -99.33,  695662,  488),
    ("Utah",           39.32, -111.09,  219882, 1860),
    ("Vermont",        44.07,  -72.67,   24906,  284),
    ("Virginia",       37.51,  -78.87,  110787,  283),
    ("Washington",     47.38, -120.45,  184661, 1011),
    ("West Virginia",  38.64,  -80.62,   62756,  472),
    ("Wisconsin",      44.27,  -89.62,  169635,  274),
    ("Wyoming",        43.00, -107.55,  253335, 2042),
]

CUTOFF_LON = -107.0
WEST_MULT  = 0.3321
EAST_MULT  = 1.0
CX_LAT, CX_LON   = 41.939, -105.082
LARAMIE_LAT, LARAMIE_LON = 41.311, -105.591

# lookup by name
STATE_MAP = {s[0]: s for s in STATES}

def eff_elev(s):
    return s[4] * (WEST_MULT if s[2] <= CUTOFF_LON else EAST_MULT)

max_eff = max(eff_elev(s) for s in STATES)
Z_SCALE = 11.0 / max_eff   # normalise so max bar ≈ 11 deg-equivalent

def state_height(name):
    s = STATE_MAP.get(name)
    return eff_elev(s) * Z_SCALE if s else 0.0

def state_color(name):
    s = STATE_MAP.get(name)
    if s is None: return 'rgba(100,100,100,0.7)'
    return ('rgba(190,55,55,0.72)' if s[2] <= CUTOFF_LON
            else 'rgba(40,100,210,0.72)')

def state_edge_color(name):
    s = STATE_MAP.get(name)
    if s is None: return 'rgba(160,160,160,0.9)'
    return ('rgba(255,140,140,0.9)' if s[2] <= CUTOFF_LON
            else 'rgba(130,180,255,0.9)')

# ── Fetch GeoJSON ──────────────────────────────────────────────────────────
print("Fetching US state boundaries...")
URL = ("https://raw.githubusercontent.com/PublicaMundi/MappingAPI"
       "/master/data/geojson/us-states.json")
with urlopen(URL) as r:
    geojson = json.load(r)
print(f"  {len(geojson['features'])} features loaded.")

# ── Triangulate a polygon ring ─────────────────────────────────────────────
def triangulate_ring(coords_2d):
    """
    Delaunay triangulation of a polygon ring.
    Returns (xs, ys, i, j, k) in deduplicated point space, or None on failure.
    Filters out large hull-bridging triangles to handle concave shapes.
    """
    if len(coords_2d) < 3:
        return None
    pts = np.unique(np.array(coords_2d, dtype=float), axis=0)
    if len(pts) < 3:
        return None
    # Check for near-collinear / degenerate set
    if np.ptp(pts[:, 0]) < 1e-6 or np.ptp(pts[:, 1]) < 1e-6:
        return None
    try:
        tri = Delaunay(pts, qhull_options='QJ')  # QJ: joggle to avoid degeneracy
    except Exception:
        return None
    simplices = tri.simplices
    areas = []
    for s in simplices:
        a, b, c = pts[s[0]], pts[s[1]], pts[s[2]]
        areas.append(abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])) / 2)
    if not areas:
        return None
    med_a = sorted(areas)[len(areas) // 2]
    threshold = med_a * 30
    good = [(s, ar) for s, ar in zip(simplices, areas) if ar < threshold and ar > 1e-12]
    if not good:
        return None
    ii = [t[0][0] for t in good]
    jj = [t[0][1] for t in good]
    kk = [t[0][2] for t in good]
    return pts[:, 0].tolist(), pts[:, 1].tolist(), ii, jj, kk

# ── Build figure ───────────────────────────────────────────────────────────
fig = go.Figure()

# Ground plane
gx = np.array([[-130, -64], [-130, -64]])
gy = np.array([[23, 23], [51, 51]])
fig.add_trace(go.Surface(
    x=gx, y=gy, z=np.zeros((2,2)),
    colorscale=[[0,'rgb(28,38,28)'],[1,'rgb(28,38,28)']],
    showscale=False, opacity=0.5, hoverinfo='skip', name='Ground'
))

# ── Per-state: top face + side walls + ground outline ─────────────────────
outline_x, outline_y = [], []          # all ground outlines in one trace
top_outline_segs = {}                  # name -> (xs, ys, zs)

SKIP_STATES = {'Alaska', 'Hawaii'}     # clip to continental US in 3D view

for feat in geojson['features']:
    name = feat['properties']['name']
    geom = feat['geometry']
    h    = state_height(name)
    face_col = state_color(name)
    edge_col = state_edge_color(name)

    rings = []
    if geom['type'] == 'Polygon':
        rings = geom['coordinates']
    elif geom['type'] == 'MultiPolygon':
        for poly in geom['coordinates']:
            rings.extend(poly)

    for ring in rings:
        coords = ring   # list of [lon, lat]
        lons = [c[0] for c in coords]
        lats = [c[1] for c in coords]
        n = len(lons)

        # ---- ground outline (all states) ----
        outline_x += lons + [None]
        outline_y += lats + [None]

        if name in SKIP_STATES:
            continue

        # ---- top-face mesh ----
        pts2d = [[c[0], c[1]] for c in coords[:-1]]   # drop repeated close point
        if len(pts2d) >= 3:
            result = triangulate_ring(pts2d)
            if result is not None:
                xs, ys, ii, jj, kk = result
                s_info = STATE_MAP.get(name)
                hover_txt = (
                    f"<b>{name}</b><br>"
                    f"Elev: {s_info[4]} m<br>"
                    f"Group: {'West x0.33' if s_info[2]<=CUTOFF_LON else 'East x1.00'}<br>"
                    f"Eff height: {eff_elev(s_info)*Z_SCALE:.2f}"
                ) if s_info else name
                fig.add_trace(go.Mesh3d(
                    x=xs, y=ys, z=[h]*len(xs),
                    i=ii, j=jj, k=kk,
                    color=face_col, flatshading=True,
                    showscale=False, showlegend=False,
                    hoverinfo='text', text=hover_txt,
                    opacity=0.88, name=name
                ))

        # ---- side walls (extrusion edges) ----
        wall_x, wall_y, wall_z = [], [], []
        step = max(1, n // 80)   # sample ~80 verticals per ring
        for idx in range(0, n, step):
            lo, la = lons[idx], lats[idx]
            wall_x += [lo, lo, None]
            wall_y += [la, la, None]
            wall_z += [0,  h,  None]
        fig.add_trace(go.Scatter3d(
            x=wall_x, y=wall_y, z=wall_z,
            mode='lines',
            line=dict(color=edge_col, width=1),
            showlegend=False, hoverinfo='skip'
        ))

        # ---- top-edge outline ----
        fig.add_trace(go.Scatter3d(
            x=lons, y=lats, z=[h]*len(lons),
            mode='lines',
            line=dict(color=edge_col, width=1.5),
            showlegend=False, hoverinfo='skip'
        ))

# ---- Ground state outlines (single trace) ----
fig.add_trace(go.Scatter3d(
    x=outline_x, y=outline_y, z=[0.02]*len(outline_x),
    mode='lines',
    line=dict(color='rgba(220,220,220,0.65)', width=1.2),
    name='State boundaries', showlegend=True, hoverinfo='skip'
))

# ── State name labels at bar tops ─────────────────────────────────────────
KEY = {"Colorado","New Mexico","Wyoming","Montana","Nevada","Utah",
       "Texas","Kansas","North Dakota","California","Idaho"}
for s in STATES:
    name, lat, lon, area, elev = s
    if name in SKIP_STATES or name not in KEY:
        continue
    h = state_height(name)
    short = name.replace("New ","N.").replace("North ","N.").replace("South ","S.")
    fig.add_trace(go.Scatter3d(
        x=[lon], y=[lat], z=[h + 0.18],
        mode='text',
        text=[f"<b>{short}</b>"],
        textfont=dict(size=9, color='white'),
        showlegend=False, hoverinfo='skip'
    ))

# ── Legend proxies ─────────────────────────────────────────────────────────
fig.add_trace(go.Scatter3d(x=[None], y=[None], z=[None], mode='markers',
    marker=dict(size=10, color='rgba(190,55,55,0.85)'),
    name='West states  (x0.33 elev)'))
fig.add_trace(go.Scatter3d(x=[None], y=[None], z=[None], mode='markers',
    marker=dict(size=10, color='rgba(40,100,210,0.85)'),
    name='East states  (x1.00 elev)'))

# ── Centroid marker ────────────────────────────────────────────────────────
cx_z_vis = max_eff * Z_SCALE * 1.55   # float well above tallest bar

fig.add_trace(go.Scatter3d(
    x=[CX_LON, CX_LON], y=[CX_LAT, CX_LAT], z=[0, cx_z_vis],
    mode='lines',
    line=dict(color='rgba(255,215,0,0.55)', width=2, dash='dash'),
    showlegend=False, hoverinfo='skip'
))
theta = np.linspace(0, 2*math.pi, 42)
fig.add_trace(go.Scatter3d(
    x=CX_LON + 0.32*np.cos(theta),
    y=CX_LAT + 0.32*np.sin(theta),
    z=np.zeros(42),
    mode='lines', line=dict(color='gold', width=2),
    showlegend=False, hoverinfo='skip'
))
fig.add_trace(go.Scatter3d(
    x=[CX_LON], y=[CX_LAT], z=[cx_z_vis],
    mode='markers+text',
    marker=dict(size=15, color='gold', symbol='diamond',
                line=dict(color='black', width=2)),
    text=['<b>Mass Centroid</b><br>41.94°N  105.08°W<br>180 km below surface'],
    textposition='top right',
    textfont=dict(size=10, color='gold'),
    name='Mass Centroid',
    hovertext='41.94N 105.08W — 180 km deep — 82 km N of Laramie',
    hoverinfo='text'
))

# ── Laramie ────────────────────────────────────────────────────────────────
fig.add_trace(go.Scatter3d(
    x=[LARAMIE_LON], y=[LARAMIE_LAT], z=[0.08],
    mode='markers+text',
    marker=dict(size=7, color='white', symbol='circle',
                line=dict(color='black', width=1.5)),
    text=['Laramie, WY'],
    textposition='top right',
    textfont=dict(size=10, color='white'),
    name='Laramie, WY',
    hovertext='Laramie WY: 41.31N 105.59W',
    hoverinfo='text'
))
fig.add_trace(go.Scatter3d(
    x=[CX_LON, LARAMIE_LON], y=[CX_LAT, LARAMIE_LAT], z=[0.08, 0.08],
    mode='lines',
    line=dict(color='rgba(255,215,0,0.45)', width=1.5, dash='dot'),
    showlegend=False, hoverinfo='skip'
))

# ── Layout ─────────────────────────────────────────────────────────────────
fig.update_layout(
    title=dict(
        text=('3D Elevation-Weighted State Map — Mass Centroid near Laramie, WY<br>'
              '<sup>Height = mean elevation &times; group multiplier &nbsp;|&nbsp; '
              '107&deg;W split &nbsp;|&nbsp; West (red) &times;0.33, '
              'East (blue) &times;1.00 &nbsp;|&nbsp; '
              'Gold diamond = centroid, 82 km N of Laramie, 180 km deep</sup>'),
        font=dict(size=14, color='white'), x=0.5
    ),
    scene=dict(
        xaxis=dict(title=dict(text='Longitude', font=dict(color='white')),
                   range=[-128, -64],
                   backgroundcolor='rgb(12,18,28)',
                   gridcolor='rgba(70,70,70,0.4)',
                   showbackground=True, tickfont=dict(color='white')),
        yaxis=dict(title=dict(text='Latitude', font=dict(color='white')),
                   range=[22, 52],
                   backgroundcolor='rgb(12,18,28)',
                   gridcolor='rgba(70,70,70,0.4)',
                   showbackground=True, tickfont=dict(color='white')),
        zaxis=dict(title=dict(text='Effective Elevation (scaled)', font=dict(color='white')),
                   backgroundcolor='rgb(8,12,22)',
                   gridcolor='rgba(70,70,70,0.4)',
                   showbackground=True, tickfont=dict(color='white')),
        bgcolor='rgb(8,12,22)',
        aspectmode='manual',
        aspectratio=dict(x=2.1, y=1.0, z=0.8),
        camera=dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0.0, y=0.0, z=-0.05),
            eye=dict(x=0.5, y=-1.85, z=1.05)
        )
    ),
    paper_bgcolor='rgb(8,12,22)',
    font=dict(color='white'),
    legend=dict(bgcolor='rgba(20,25,40,0.85)',
                bordercolor='rgba(120,120,120,0.4)', borderwidth=1,
                font=dict(size=11, color='white')),
    margin=dict(l=0, r=0, t=75, b=0),
    width=1440, height=880
)

out = 'centroid_3d.html'
fig.write_html(out, include_plotlyjs='cdn')
print(f"Saved: {out}")
