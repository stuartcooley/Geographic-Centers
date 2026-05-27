#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
US Geographic Centroid - East/West elevation-thickness split
Goal: find longitude cutoff + west/east multiplier ratio that places
the area*elevation-weighted centroid near Laramie, WY (41.311N, 105.591W)
"""
import sys, math

# Force UTF-8 output on Windows
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np

# (name, lat_centroid, lon_centroid, land_area_km2, mean_elevation_m)
# Areas: US Census Bureau 2020; Centroids: USGS; Elevations: USGS/NOAA
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

TARGET_LAT, TARGET_LON = 41.311, -105.591
R = 6371.0

def unit_vec(lat_d, lon_d):
    la, lo = math.radians(lat_d), math.radians(lon_d)
    return np.array([math.cos(la)*math.cos(lo),
                     math.cos(la)*math.sin(lo),
                     math.sin(la)])

def centroid(wt_list):
    C = np.zeros(3)
    total = 0.0
    for lat, lon, w in wt_list:
        C += w * unit_vec(lat, lon)
        total += w
    C /= total
    mag = float(np.linalg.norm(C))
    n = C / mag
    return (float(np.degrees(np.arcsin(n[2]))),
            float(np.degrees(np.arctan2(n[1], n[0]))),
            float(R * (1 - mag)), mag)

def hav(la1, lo1, la2, lo2):
    a = (math.sin(math.radians(la2-la1)/2)**2 +
         math.cos(math.radians(la1))*math.cos(math.radians(la2)) *
         math.sin(math.radians(lo2-lo1)/2)**2)
    return 2*R*math.asin(math.sqrt(a))

def run(cut, wmult, emult=1.0):
    wts = [(s[1], s[2], s[3]*s[4]*(wmult if s[2] <= cut else emult))
           for s in STATES]
    return centroid(wts)

# ---------- baselines ----------
wts0 = [(s[1], s[2], float(s[3])) for s in STATES]
lc0, lo0, d0, m0 = centroid(wts0)
print(f"Baseline area-only:                 {lc0:.3f}N {abs(lo0):.3f}W  depth {d0:.0f} km")

wts1 = [(s[1], s[2], float(s[3]*s[4])) for s in STATES]
lc1, lo1, d1, m1 = centroid(wts1)
print(f"Elevation-weighted (no split):      {lc1:.3f}N {abs(lo1):.3f}W  depth {d1:.0f} km")
print(f"Target Laramie WY:                  {TARGET_LAT:.3f}N {abs(TARGET_LON):.3f}W")
print()

# ---------- Grid search ----------
print("Searching cutoffs 107W to 89W, ratios 0.01-100 (log scale)...")
ratios = np.exp(np.linspace(math.log(0.01), math.log(100), 1200))
best = {"dist": 1e9}
results_table = []

for cut_deg in range(107, 88, -1):
    cut = -float(cut_deg)
    row_best = {"dist": 1e9}
    for r in ratios:
        lc, loc, dep, mag = run(cut, float(r))
        d = hav(lc, loc, TARGET_LAT, TARGET_LON)
        if d < row_best["dist"]:
            row_best = {"dist": d, "r": float(r), "lat": lc, "lon": loc,
                        "dep": dep, "mag": mag, "cut": cut_deg}
    results_table.append(row_best)
    if row_best["dist"] < best["dist"]:
        best = row_best

# ---------- Progression table ----------
print(f"\n{'Cutoff W':>9} {'W/E ratio':>10} {'Lat N':>7} {'Lon W':>7} {'Depth km':>9} {'Dist km':>9}")
print("-" * 58)
for row in results_table:
    marker = " <<< BEST" if row is best else ""
    print(f"  {row['cut']:>5}W   {row['r']:>10.4f} {row['lat']:>7.3f} {abs(row['lon']):>7.3f}"
          f" {row['dep']:>9.1f} {row['dist']:>9.1f}{marker}")

# ---------- Best solution detail ----------
print(f"\n{'='*65}")
print(f"OPTIMAL SOLUTION")
print(f"  Longitude cutoff  : {best['cut']}W")
west = [s[0] for s in STATES if s[2] <= -best['cut']]
east = [s[0] for s in STATES if s[2] >  -best['cut']]
print(f"  Western states ({len(west)}): {', '.join(west)}")
print(f"  Eastern states ({len(east)}): {', '.join(east)}")
print(f"  W/E multiplier    : {best['r']:.4f}  "
      f"(west elev weight x{best['r']:.4f}, east x1.0)")
print(f"  Centroid          : {best['lat']:.3f}N, {abs(best['lon']):.3f}W")
print(f"  Distance to Laramie: {best['dist']:.1f} km")
print()
print(f"  Depth below surface:")
print(f"    |C_bar| = {best['mag']:.6f}")
print(f"    d = R x (1 - |C_bar|) = {R:.0f} x (1 - {best['mag']:.6f})"
      f" = {R*(1-best['mag']):.1f} km")

# ---------- Per-state weight table at optimal ----------
print(f"\nState weights at optimal {best['cut']}W cutoff  (W/E = {best['r']:.4f}):")
print(f"  {'State':<20} {'Grp':>4} {'Elev':>6} {'Area km2':>10} {'Eff weight':>14}")
print("  " + "-"*58)
cut_opt = -float(best['cut'])
r_opt = best['r']
for s in sorted(STATES, key=lambda x: x[2]):
    grp = "WEST" if s[2] <= cut_opt else "east"
    mult = r_opt if s[2] <= cut_opt else 1.0
    ew = s[3]*s[4]*mult
    print(f"  {s[0]:<20} {grp:>4} {s[4]:>6} {s[3]:>10,} {ew:>14,.0f}")
