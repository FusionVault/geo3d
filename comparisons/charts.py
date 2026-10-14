#!/usr/bin/env python3
"""Draw the comparison charts (SVG, light and dark) from the suite's own results.

    cargo run --release --bin accuracy > results/accuracy.md
    cargo bench                       # writes target/criterion/
    python3 charts.py                 # → results/speed.csv, charts/*.svg and the READMEs' alt texts

Standard library only. The charts use the emphasis form: geo3d in the accent colour, every other
crate in grey, each bar labelled, so identity never depends on colour alone. Errors that span many
orders of magnitude use a log scale. The tables behind every chart are in README.md.
"""

import csv
import json
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
CHARTS = os.path.join(HERE, "charts")
A_WGS84 = 6_378_137.0

THEMES = {
    "light": dict(surface="#fcfcfb", border="#e1e0d9", ink="#0b0b0b", ink2="#52514e", muted="#6f6d68",
                  grid="#e9e8e2", axis="#c3c2b7", accent="#2a78d6", other="#898781",
                  good="#0ca30c", bad="#d03b3b"),
    "dark": dict(surface="#1a1a19", border="#2c2c2a", ink="#ffffff", ink2="#c3c2b7", muted="#9b9a93",
                 grid="#2c2c2a", axis="#383835", accent="#3987e5", other="#75746d",
                 good="#0ca30c", bad="#d03b3b"),
}
FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"
CRATES = ["geo3d", "map_3d", "nav-types", "sguaba", "jord", "geographiclib-rs", "geo"]


# ── Reading the results ────────────────────────────────────────────────────────────────────────

def tables(markdown):
    """`{heading: [row dicts]}` for every Markdown table under a `###` heading."""
    out, heading, header = {}, None, None
    for line in markdown.splitlines():
        if line.startswith("### "):
            heading, header = line[4:].strip(), None
        elif line.startswith("|") and heading:
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if header is None:
                header = cells
            elif not set(line.replace("|", "").strip()) <= set("-: "):
                out.setdefault(heading, []).append(dict(zip(header, cells)))
    return out


def find(tabs, prefix):
    return next(rows for name, rows in tabs.items() if name.startswith(prefix))


def num(cell):
    m = re.search(r"-?[\d.]+(?:e-?\d+)?", cell.replace(" ", ""))
    return float(m.group()) if m else None


def speed_rows():
    """ns per call by operation and crate, from criterion (per 1 000 calls) or results/speed.csv."""
    root = os.path.join(HERE, "target", "criterion")
    path = os.path.join(RESULTS, "speed.csv")
    if os.path.isdir(root):
        rows = []
        for op in sorted(os.listdir(root)):
            for crate in CRATES:
                est = os.path.join(root, op, crate, "new", "estimates.json")
                if os.path.exists(est):
                    with open(est) as f:
                        rows.append((op, crate, json.load(f)["median"]["point_estimate"] / 1000.0))
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["operation", "crate", "ns_per_call"])
            w.writerows((op, c, f"{ns:.1f}") for op, c, ns in rows)
    with open(path) as f:
        return [(r["operation"], r["crate"], float(r["ns_per_call"])) for r in csv.DictReader(f)]


# ── Formatting ─────────────────────────────────────────────────────────────────────────────────

def sig(v):
    return f"{v:.0f}" if v >= 10 else f"{v:.1f}" if v >= 1 else f"{v:.2f}"


def metres(x):
    for scale, unit in [(1e3, "km"), (1.0, "m"), (1e-2, "cm"), (1e-3, "mm"), (1e-6, "µm"), (1e-9, "nm")]:
        if x >= scale or unit == "nm":
            return f"{sig(x / scale)} {unit}"


def mix(fg, bg, alpha):
    """`fg` at `alpha` over `bg`, as a solid colour (some renderers ignore SVG opacity)."""
    f, b = (tuple(int(c[i:i + 2], 16) for i in (1, 3, 5)) for c in (fg, bg))
    return "#" + "".join(f"{round(x * alpha + y * (1 - alpha)):02x}" for x, y in zip(f, b))


def decade_label(e):
    units = {-9: "1 nm", -6: "1 µm", -3: "1 mm", 0: "1 m", 3: "1 km", 6: "1000 km"}
    if e in units:
        return units[e]
    base = 3 * math.floor(e / 3)
    return units.get(base, "").replace("1 ", f"{10 ** (e - base)} ")


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text_w(s, size):
    return len(s) * size * 0.56  # a fair estimate for a system sans


# ── SVG pieces ─────────────────────────────────────────────────────────────────────────────────

class Svg:
    def __init__(self, width, theme):
        self.w, self.t, self.parts, self.h = width, THEMES[theme], [], 0

    def add(self, s):
        self.parts.append(s)

    def text(self, x, y, s, size=12.5, fill="ink2", weight=400, anchor="start"):
        self.add(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" font-weight="{weight}" '
                 f'fill="{self.t[fill]}" text-anchor="{anchor}">{esc(s)}</text>')

    def header(self, title, subtitle):
        self.text(24, 34, title, 16, "ink", 600)
        y = 56
        for line in subtitle:
            self.text(24, y, line, 12.5, "ink2")
            y += 18
        return y + 6

    def bar(self, x0, y, x1, h, fill):
        """A horizontal bar: square at the baseline, 4 px rounded at the data end."""
        r = min(4.0, (x1 - x0) / 2, h / 2)
        self.add(f'<path d="M{x0:.1f},{y:.1f} H{x1 - r:.1f} A{r},{r} 0 0 1 {x1:.1f},{y + r:.1f} '
                 f'V{y + h - r:.1f} A{r},{r} 0 0 1 {x1 - r:.1f},{y + h:.1f} H{x0:.1f} Z" fill="{self.t[fill]}"/>')

    def status(self, cx, cy, ok):
        """A status icon: a filled circle with a check or a cross (never colour alone)."""
        self.add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="7" fill="{self.t["good" if ok else "bad"]}"/>')
        if ok:
            d = f"M{cx - 3.2:.1f},{cy + 0.2:.1f} L{cx - 0.9:.1f},{cy + 2.6:.1f} L{cx + 3.4:.1f},{cy - 2.6:.1f}"
        else:
            d = (f"M{cx - 2.6:.1f},{cy - 2.6:.1f} L{cx + 2.6:.1f},{cy + 2.6:.1f} "
                 f"M{cx + 2.6:.1f},{cy - 2.6:.1f} L{cx - 2.6:.1f},{cy + 2.6:.1f}")
        self.add(f'<path d="{d}" stroke="#ffffff" stroke-width="1.8" fill="none" stroke-linecap="round" '
                 f'stroke-linejoin="round"/>')

    def save(self, name, theme, alt):
        t = self.t
        body = "\n".join(self.parts)
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" height="{self.h}" '
               f'viewBox="0 0 {self.w} {self.h}" font-family="{FONT}" role="img" aria-label="{esc(alt)}">\n'
               f'<title>{esc(alt)}</title>\n'
               f'<rect x="0.5" y="0.5" width="{self.w - 1}" height="{self.h - 1}" rx="10" fill="{t["surface"]}" '
               f'stroke="{t["border"]}"/>\n{body}\n</svg>\n')
        with open(os.path.join(CHARTS, f"{name}-{theme}.svg"), "w") as f:
            f.write(svg)


def bars_log(name, title, subtitle, rows, footnote, alt, width=760):
    """Horizontal bars on a log axis in metres; rows are (label, detail, value, is_geo3d)."""
    for theme in THEMES:
        s = Svg(width, theme)
        y = s.header(title, subtitle)
        label_w = max(text_w(lbl, 13) for lbl, _, _, _ in rows) + 8
        detail_w = max(text_w(d, 11.5) for _, d, _, _ in rows)
        x0 = 24 + max(label_w, detail_w) + 12
        x1 = width - 96
        lo = math.floor(math.log10(min(v for *_, v, _ in rows)))
        hi = math.ceil(math.log10(max(v for *_, v, _ in rows)))
        if hi - lo > 6:  # snap a long axis to unit boundaries: 1 nm, 1 µm, 1 mm, 1 m, 1 km
            lo, hi = 3 * math.floor(lo / 3), 3 * math.ceil(hi / 3)
        scale = lambda v: x0 + (math.log10(v) - lo) / (hi - lo) * (x1 - x0)
        row_h, bar_h = 40, 18
        plot_top, plot_bot = y + 8, y + 8 + row_h * len(rows)
        step = 3 if hi - lo > 6 else 1
        for e in range(lo, hi + 1):
            if (e - lo) % step == 0:
                x = scale(10 ** e)
                s.add(f'<line x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" y2="{plot_bot}" stroke="{s.t["grid"]}"/>')
                s.text(x, plot_bot + 16, decade_label(e), 11, "muted", anchor="middle")
        s.add(f'<line x1="{x0}" y1="{plot_top}" x2="{x0}" y2="{plot_bot}" stroke="{s.t["axis"]}"/>')
        for i, (lbl, detail, v, ours) in enumerate(rows):
            cy = plot_top + i * row_h + row_h / 2
            s.text(24, cy - 1, lbl, 13, "ink", 600 if ours else 400)
            s.text(24, cy + 13, detail, 11.5, "muted")
            s.bar(x0, cy - bar_h / 2, max(scale(v), x0 + 3), bar_h, "accent" if ours else "other")
            s.text(max(scale(v), x0 + 3) + 7, cy + 4.5, metres(v), 12.5, "ink", 600 if ours else 400)
        y = plot_bot + 40
        for line in footnote:
            s.text(24, y, line, 11.5, "muted")
            y += 16
        s.h = int(y + 6)
        s.save(name, theme, alt)


def scorecard(name, title, subtitle, columns, rows, footnote, alt, kind, width=900):
    """A grid of crates (columns) by cases (rows); kind "status" draws ✓/✗ icons with labels."""
    for theme in THEMES:
        s = Svg(width, theme)
        y = s.header(title, subtitle)
        label_w = max(text_w(r[0], 12.5) for r in rows) + 16
        col_w = (width - 24 - 24 - label_w) / len(columns)
        cx0 = 24 + label_w
        top = y + 4
        # geo3d's column carries the emphasis wash.
        wash = mix(s.t["accent"], s.t["surface"], 0.10 if theme == "light" else 0.18)
        s.add(f'<rect x="{cx0:.1f}" y="{top:.1f}" width="{col_w:.1f}" height="{30 + 30 * len(rows):.1f}" '
              f'rx="6" fill="{wash}"/>')
        for j, c in enumerate(columns):
            s.text(cx0 + col_w * (j + 0.5), top + 20, c, 12, "ink", 600 if c == "geo3d" else 500, "middle")
        for i, (label, cells) in enumerate(rows):
            ry = top + 30 + 30 * i
            s.add(f'<line x1="24" y1="{ry:.1f}" x2="{width - 24}" y2="{ry:.1f}" stroke="{s.t["grid"]}"/>')
            s.text(24, ry + 19, label, 12.5, "ink")
            for j, cell in enumerate(cells):
                cx = cx0 + col_w * (j + 0.5)
                if cell is None:
                    s.text(cx, ry + 19, "—", 12.5, "muted", anchor="middle")
                elif kind == "status":
                    ok, word = cell
                    w = 18 + text_w(word, 11.5)
                    s.status(cx - w / 2 + 7, ry + 15, ok)
                    s.text(cx - w / 2 + 18, ry + 19, word, 11.5, "ink2")
                elif cell == "✓":
                    # Drawn, not typed: not every font has a check mark glyph.
                    s.add(f'<path d="M{cx - 5:.1f},{ry + 15:.1f} L{cx - 1.5:.1f},{ry + 18.5:.1f} '
                          f'L{cx + 5.5:.1f},{ry + 10.5:.1f}" stroke="{s.t["ink"]}" stroke-width="2" fill="none" '
                          f'stroke-linecap="round" stroke-linejoin="round"/>')
                else:
                    s.text(cx, ry + 19, cell, 12, "ink2", 500, "middle")
        y = top + 30 + 30 * len(rows) + 22
        for line in footnote:
            s.text(24, y, line, 11.5, "muted")
            y += 16
        s.h = int(y + 6)
        s.save(name, theme, alt)


def speed_panels(name, title, subtitle, panels, footnote, alt, width=860):
    """Small multiples: one panel of horizontal bars (ns per call) per operation."""
    cols, gap = 2, 28
    pw = (width - 48 - gap) / cols
    for theme in THEMES:
        s = Svg(width, theme)
        y = s.header(title, subtitle)
        row_h, bar_h = 22, 13
        heights = [34 + row_h * len(p[1]) for p in panels]
        for k, (ptitle, rows) in enumerate(panels):
            col, rowi = k % cols, k // cols
            px = 24 + col * (pw + gap)
            py = y + sum(max(heights[r * cols:(r + 1) * cols]) + 18 for r in range(rowi))
            s.text(px, py + 14, ptitle, 13, "ink", 600)
            name_w = 96
            bx0, bx1 = px + name_w, px + pw - 52
            vmax = max(v for _, v, _ in rows)
            for i, (crate, v, ours) in enumerate(rows):
                cy = py + 30 + i * row_h + row_h / 2
                s.text(px, cy + 4, crate, 11.5, "ink" if ours else "ink2", 600 if ours else 400)
                x = bx0 + v / vmax * (bx1 - bx0)
                s.bar(bx0, cy - bar_h / 2, x, bar_h, "accent" if ours else "other")
                s.text(x + 6, cy + 4, f"{v / 1000:.2f} µs" if v >= 1000 else f"{v:.0f} ns", 11.5, "ink",
                       600 if ours else 400)
            s.add(f'<line x1="{bx0}" y1="{py + 30}" x2="{bx0}" y2="{py + 30 + row_h * len(rows)}" '
                  f'stroke="{s.t["axis"]}"/>')
        rows_n = math.ceil(len(panels) / cols)
        y += sum(max(heights[r * cols:(r + 1) * cols]) + 18 for r in range(rows_n)) + 8
        for line in footnote:
            s.text(24, y, line, 11.5, "muted")
            y += 16
        s.h = int(y + 6)
        s.save(name, theme, alt)


# ── The charts ─────────────────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(CHARTS, exist_ok=True)
    with open(os.path.join(RESULTS, "accuracy.md")) as f:
        tabs = tables(f.read())

    # 1. Round trip: worst position error over every height, as metres along the surface/orbit.
    lat, lon, hgt = (find(tabs, p) for p in ("Max latitude", "Max longitude", "Max height"))
    worst = {}
    for rl, ro, rh in zip(lat, lon, hgt):
        h = num(rl["Height"]) * 1000.0
        radius = A_WGS84 + h
        for c in ("geo3d", "map_3d", "nav-types", "sguaba", "jord"):
            err = max(math.radians(num(rl[c])) * radius, math.radians(num(ro[c])) * radius, num(rh[c]))
            worst[c] = max(worst.get(c, 0.0), err)
    bars_log(
        "roundtrip", "Round trip: latitude/longitude/height → x, y, z → back",
        ["Worst position error over 66 960 points per height, from 11 km underground to geostationary",
         "orbit. Each crate uses its own conversions. Log scale; shorter is better."],
        [("geo3d", "closed form, then one correction step", worst["geo3d"], True),
         ("nav-types", "closed form", worst["nav-types"], False),
         ("jord", "closed form", worst["jord"], False),
         ("sguaba", "closed form", worst["sguaba"], False),
         ("map_3d", "one fixed correction step", worst["map_3d"], False)],
        ["map_3d's error grows with height: about 30 m at GPS altitude and 80 m at geostationary orbit."],
        f"Worst round-trip position error up to geostationary height: geo3d {metres(worst['geo3d'])}, "
        f"map_3d {metres(worst['map_3d'])}")

    # 2. Surface distance: the best each crate offers, against GeographicLib.
    inv = find(tabs, "Inverse")
    geo3d_inv = max(num(r["Max distance error (m)"]) for r in inv if r["Library"] == "geo3d")
    gc = {r["Library"]: num(r["Max error (km)"]) * 1000.0 for r in find(tabs, "Great-circle")}
    bars_log(
        "surface-distance", "Distance along the surface: the best each crate offers",
        ["Worst error against GeographicLib over 110 000 point pairs, 10 000 of them nearly",
         "antipodal. Log scale; shorter is better."],
        [("geo3d", "Karney geodesic on the ellipsoid", geo3d_inv, True),
         ("map_3d", "great circle on a sphere", gc["map_3d"], False),
         ("jord", "great circle on a sphere", gc["jord"], False),
         ("nav-types", "great circle on a sphere", gc["nav-types"], False),
         ("sguaba", "great circle on a sphere", gc["sguaba"], False)],
        ["geographiclib-rs and geo's Geodesic are GeographicLib itself, the reference: geo3d agrees with them",
         "to 4 nm and, like them, answers every pair. geo3d also has a great-circle helper for quick estimates."],
        f"Worst surface-distance error: geo3d {metres(geo3d_inv)}, great-circle crates "
        f"{metres(min(gc.values()))} to {metres(max(gc.values()))}")

    # 3. Look angles on the hard cases: straight up and on/near the cardinal directions.
    hard = {}
    for prefix in ("2400 pairs", "2000 pairs"):
        for r in find(tabs, prefix):
            v = num(r["Look angles vs the built pointing (µm)"].replace(" ", "")) if r[
                "Look angles vs the built pointing (µm)"] != "—" else None
            if v is not None:
                hard[r["Library"]] = max(hard.get(r["Library"], 0.0), v * 1e-6)
    bars_log(
        "look-angles", "Look angles on the hard cases",
        ["Worst error for targets straight up and just off it, and on or just off due north, east, south",
         "and west, out to geostationary range. Log scale; shorter is better."],
        [("geo3d", "atan2 for both angles", hard["geo3d"], True),
         ("map_3d", "atan2 for both angles", hard["map_3d"], False),
         ("jord", "asin for elevation", hard["jord"], False),
         ("sguaba", "acos for both angles", hard["sguaba"], False)],
        ["nav-types has no look angles. Error is how far the answer moves the target."],
        f"Worst look-angle error on hard cases: geo3d {metres(hard['geo3d'])}, jord {metres(hard['jord'])}, "
        f"sguaba {metres(hard['sguaba'])}")

    # 4. Robustness: what each crate returns on the edge cases.
    rob = find(tabs, "What each library returns")
    short = {
        "ECEF → geodetic at the Earth's centre": "At the Earth's centre",
        "ECEF → geodetic exactly above the North Pole": "Exactly over the North Pole",
        "Round trip exactly at the South Pole": "Exactly at the South Pole",
        "ENU across the antimeridian": "Across the antimeridian",
        "Round trip at 1e8 m altitude": "At 100 000 km altitude",
        "Geodetic → ECEF with a NaN latitude": "NaN latitude in",
        "ECEF → geodetic with a NaN x": "NaN x, y, z in",
        "Geodesic between exact antipodes": "Geodesic, exact antipodes",
        "Geodesic with a NaN latitude": "Geodesic, NaN latitude in",
    }

    def classify(case, cell):
        if cell == "—":
            return None
        nan_in = "NaN" in case
        if cell.startswith("ok"):
            return (True, "correct")
        if nan_in and (cell == "NaN" or cell.startswith("rejected")):
            return (True, "NaN out" if cell == "NaN" else "rejected")
        if "panic" in cell:
            return (False, "panic")
        if cell.startswith("wrong"):
            return (False, "wrong")
        return (False, "NaN")

    rows = []
    for r in rob:
        case = r["Case"]
        label = next(v for k, v in short.items() if case.startswith(k))
        rows.append((label, [classify(case, r[c]) for c in CRATES]))
    scorecard(
        "robustness", "Edge cases: what each crate returns",
        ["A correct answer, or for NaN input NaN out or a clean rejection, counts as a pass. A panic, a",
         "wrong answer, or NaN from valid input is a failure."],
        CRATES, rows, ["— : the crate does not offer the operation."],
        "Edge-case scorecard: geo3d passes " + str(sum(1 for _, c in rows if c[0] and c[0][0])) + " of "
        + str(len(rows)) + " cases", "status")

    # 5. Polygon areas: a scorecard, since geo3d and geo's GeodesicArea match the reference exactly.
    area = find(tabs, "Area of a polygon")
    area_cols = [("geo3d", "geo3d"), ("geo (`GeodesicArea`)", "geo, geodesic"),
                 ("geo (`ChamberlainDuquetteArea`, sphere)", "geo, sphere"), ("jord (spherical excess, sphere)", "jord, sphere")]
    area_rows = {"fields": "Fields, 10 m to 1 km", "regions": "Regions, 1 km to 100 km",
                 "across": "Across the antimeridian", "around": "Around a pole"}

    def area_cell(r):
        m2, share = num(r["Max error (m²)"]), r["Max error (share of the area)"]
        if int(r["Failures"]) > 0:
            return (False, "fails")
        if m2 <= 0.1:
            return (True, "matches")
        pc = num(share.replace("%", "").strip())
        return (False, "wrong" if pc > 10 else f"{pc:.1f} % off")

    rows = []
    for key, label in area_rows.items():
        cells = [area_cell(next(r for r in area if r["Library"] == lib and r["Polygons"].startswith(key)))
                 for lib, _ in area_cols]
        rows.append((label, cells))
    worst_g3 = max(num(r["Max error (m²)"]) for r in area if r["Library"] == "geo3d")
    scorecard(
        "polygon-area", "Polygon area: does it match GeographicLib?",
        ["1 000 polygons of each kind, 3 to 32 vertices. A pass is within 0.1 m² of GeographicLib, the",
         "reference. Spherical methods treat the Earth as a ball, so they are 0.7 % to 1.4 % off."],
        [c for _, c in area_cols], rows,
        [f"geo3d's largest difference from GeographicLib is {worst_g3:g} m², around a pole. geo's geodesic area calls",
         "geographiclib-rs. \"wrong\": more than 10 % off; geo's sphere method can't cross ±180° or circle a pole."],
        "Polygon area scorecard: geo3d matches GeographicLib on all four kinds of polygon; spherical methods are "
        "0.7 % to 1.4 % off, and geo's fails across the antimeridian and around the poles", "status", width=820)

    # 6. Cross-track distance on short paths, against the exact answer.
    xt = {r["Library"]: num(r["Max error (m)"]) for r in find(tabs, "Distance from a point to a path")
          if r["Paths"].startswith("short")}
    textbook = next(k for k in xt if k.startswith("textbook"))
    bars_log(
        "cross-track", "How far off the path? Short paths",
        ["Worst error over 10 000 paths from 1 µm to 1 km long, with a point 1 m to 1 km to the side,",
         "against the exact answer. Log scale; shorter is better."],
        [("geo3d", "n-vectors, from the latitude/longitude differences", xt["geo3d"], True),
         ("geo", "haversine and bearings", xt["geo"], False),
         ("jord", "n-vectors, S2's stable cross product", xt["jord"], False),
         ("textbook", "n-vectors, n1 × n2 (NavLab, nvector)", xt[textbook], False)],
        ["On paths of 1 km to 1 000 km, geo3d, geo and jord are all within a few nanometres; the textbook",
         "formula is within 2 µm. It loses precision as the path gets shorter."],
        f"Worst cross-track error on short paths: geo3d {metres(xt['geo3d'])}, geo {metres(xt['geo'])}, "
        f"jord {metres(xt['jord'])}, textbook n-vector formula {metres(xt[textbook])}")

    # 7. Attitude near straight up: how far the rebuilt orientation points off, at 1 000 km.
    att = find(tabs, "Yaw, pitch and roll")
    near = lambda lib, prefix: max(num(r["Max turn error at 1 000 km (m)"]) for r in att
                                   if r["Library"] == lib and r["Attitudes"].startswith(prefix))
    g3 = max(near("geo3d", "within"), near("geo3d", "via ECEF"))
    away = max(num(r["Max turn error at 1 000 km (m)"]) for r in att if r["Attitudes"] == "any pitch")
    bars_log(
        "attitude", "Yaw, pitch and roll near straight up",
        ["Worst error over 10 000 attitudes within 1e-15° to 0.1° of ±90° pitch, as how far the rebuilt",
         "orientation points off at 1 000 km. Log scale; shorter is better."],
        [("geo3d", "its own round trip, and via ECEF", g3, True),
         ("sguaba", "its own round trip (asin, via nalgebra)", near("sguaba", "within"), False),
         ("jord", "read back after a trip through ECEF", near("jord", "via ECEF"), False)],
        [f"Away from the vertical, all three are within {metres(away)} at 1 000 km. jord's own round trip is",
         "exact near the vertical too; it fails when the matrix carries ordinary rounding, as one from ECEF does."],
        f"Worst attitude error near ±90° pitch, at 1 000 km: geo3d {metres(g3)}, sguaba "
        f"{metres(near('sguaba', 'within'))}, jord via ECEF {metres(near('jord', 'via ECEF'))}")

    # 8. Features (from comparisons/README.md, "Features and footprint").
    Y, N = "✓", None
    feats = [
        ("No dependencies", [Y, Y, N, N, Y, N, N]),
        ("No unsafe code (forbidden)", [Y, N, N, N, Y, N, N]),
        ("Runs without std (no_std)", [Y, N, N, Y, N, N, N]),
        ("Ellipsoidal geodesics", [Y, N, N, N, N, Y, Y]),
        ("Geodesic polygon area", [Y, N, N, N, "sphere", Y, Y]),
        ("Uncertainty (covariance)", [Y, N, N, N, N, N, N]),
        ("Look angles", [Y, Y, N, Y, Y, N, N]),
        ("Ray–ellipsoid ground strike", [Y, N, N, N, N, N, N]),
        ("3D volumes", [Y, N, N, N, N, N, N]),
        ("3D closest approach", [Y, N, N, N, "sphere", N, N]),
        ("Datum shift (Helmert)", [Y, N, N, N, N, N, N]),
        ("Inertial frame (TEME)", [Y, Y, N, N, N, N, N]),
        ("Any ellipsoid", [Y, "10 named", N, N, Y, Y, N]),
        ("Path tools (cross-track, crossings, centre)", [Y, N, N, N, Y, N, "some"]),
        ("Attitude (yaw, pitch, roll)", [Y, N, N, Y, Y, N, N]),
    ]
    scorecard(
        "features", "Features",
        ["What each crate offers, checked in its source at the versions compared."],
        CRATES, feats, ["— : not offered."],
        "Feature matrix: geo3d is the only crate with ground strike, 3D volumes, closest approach, covariance and "
        "datum shifts",
        "text")

    # 9. Speed: small multiples, one panel per operation, fixed crate order.
    sp = speed_rows()
    ops = [
        ("geodetic_to_ecef", "Latitude/longitude/height → x, y, z"),
        ("ecef_to_geodetic", "x, y, z → latitude/longitude/height"),
        ("enu_one_observer", "ENU, frame reused"),
        ("enu_single_call", "ENU, single call"),
        ("aer_one_observer", "Look angles, frame reused"),
        ("aer_single_call", "Look angles, single call"),
        ("geodesic_direct", "Geodesic direct"),
        ("geodesic_inverse", "Geodesic inverse"),
        ("great_circle_distance", "Great-circle distance"),
        ("ecef_to_geodetic_orbit", "x, y, z → lat/lon/height, orbit heights"),
        ("polygon_area", "Polygon area, 10 vertices"),
    ]
    panels = []
    for op, ptitle in ops:
        rows = [(c + ("*" if c == "geo" and op.startswith("geodesic") else ""), ns, c == "geo3d")
                for c in CRATES for o, cc, ns in sp if o == op and cc == c]
        if rows:
            panels.append((ptitle, rows))
    # Within 3 % of the fastest is a tie: closer than run-to-run noise.
    best = [(min(v for _, v, _ in rows), next(v for _, v, ours in rows if ours)) for _, rows in panels]
    def runner_up(rows):
        return sorted(v for _, v, _ in rows)[1] if len(rows) > 1 else float("inf")
    outright = sum(1 for (lo, ours), (_, rows) in zip(best, panels) if ours == lo and runner_up(rows) > ours * 1.03)
    tied = sum(1 for (lo, ours), (_, rows) in zip(best, panels)
               if (ours == lo and runner_up(rows) <= ours * 1.03) or lo < ours <= lo * 1.03)
    speed_panels(
        "speed", "Time per call",
        ["Median nanoseconds per call on an AMD Ryzen 9 5950X, rustc 1.98, no target-cpu=native. Shorter",
         "is faster. Each panel has its own scale."],
        panels,
        ["* geo's Geodesic returns only the distance (inverse) or the position (direct); geo3d and geographiclib-rs",
         "also return the azimuths."],
        f"Time per call for {len(panels)} operations; geo3d is fastest in {outright} and tied in {tied}")


def sync_alt_texts():
    """Copy each chart's alt text into the READMEs that show it, so the two never drift apart."""
    for path, prefix in ((os.path.join(HERE, "..", "README.md"), "comparisons/charts/"),
                         (os.path.join(HERE, "README.md"), "charts/")):
        with open(path) as f:
            text = f.read()

        def alt(m):
            with open(os.path.join(CHARTS, m.group(1) + "-light.svg")) as f:
                label = re.search(r'aria-label="([^"]*)"', f.read()).group(1)
            return f'<img alt="{label}" src="{prefix}{m.group(1)}-light.svg">'

        with open(path, "w") as f:
            f.write(re.sub(r'<img alt="[^"]*" src="' + re.escape(prefix) + r'([a-z-]+)-light\.svg">', alt, text))


if __name__ == "__main__":
    main()
    sync_alt_texts()
