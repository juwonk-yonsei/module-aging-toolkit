"""Five-way overlap check for the v15 figures.

Text-versus-text uses the ink polygon, not the axis-aligned box of a rotated
label. The other four checks use window extents: text crossing into another
axes, text whose anchor is inside an axes but whose ink crosses that spine,
text or legend boxes covering data, and text outside the canvas.
"""
from __future__ import annotations

import numpy as np
from matplotlib.transforms import Bbox


def _ink(text, renderer):
    from matplotlib.text import Text
    # Annotation.get_window_extent also covers the arrow. The letters are the Text box.
    ext = Text.get_window_extent(text, renderer)
    poly = np.array([[ext.x0, ext.y0], [ext.x1, ext.y0], [ext.x1, ext.y1], [ext.x0, ext.y1]])
    return ext, poly


def _poly_overlap(a, b, min_area=2.0):
    if a is None or b is None:
        return False
    a_min, a_max = a.min(axis=0), a.max(axis=0)
    b_min, b_max = b.min(axis=0), b.max(axis=0)
    if a_max[0] <= b_min[0] or b_max[0] <= a_min[0] or a_max[1] <= b_min[1] or b_max[1] <= a_min[1]:
        return False
    area = (min(a_max[0], b_max[0]) - max(a_min[0], b_min[0])) * (min(a_max[1], b_max[1]) - max(a_min[1], b_min[1]))
    if area < min_area:
        return False
    axes = np.stack([a[1] - a[0], a[2] - a[1], b[1] - b[0], b[2] - b[1]])
    for axis in axes:
        norm = np.linalg.norm(axis)
        if norm == 0:
            continue
        axis = axis / norm
        pa = a @ axis
        pb = b @ axis
        if pa.max() <= pb.min() + 0.4 or pb.max() <= pa.min() + 0.4:
            return False
    return True


def _point_in_poly(pt, poly):
    x, y = pt
    inside = False
    j = len(poly) - 1
    for i in range(len(poly)):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / ((yj - yi) or 1e-9) + xi):
            inside = not inside
        j = i
    return inside


def _segments_cross(a, b, c, d):
    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    return orient(a, b, c) * orient(a, b, d) < 0 and orient(c, d, a) * orient(c, d, b) < 0


def _segment_hits_poly(p0, p1, poly):
    if _point_in_poly(p0, poly) or _point_in_poly(p1, poly):
        return True
    for i in range(len(poly)):
        if _segments_cross(p0, p1, poly[i], poly[(i + 1) % len(poly)]):
            return True
    return False


def _bb_area(a, b):
    if not a.overlaps(b):
        return 0.0
    x0, x1 = max(a.x0, b.x0), min(a.x1, b.x1)
    y0, y1 = max(a.y0, b.y0), min(a.y1, b.y1)
    return max(0.0, x1 - x0) * max(0.0, y1 - y0)


def _collect(fig, renderer):
    items = []
    for ax in fig.axes:
        artists = list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]
        if ax.axison:
            artists += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        legend = ax.get_legend()
        if legend is not None:
            artists += list(legend.get_texts())
        for t in artists:
            label = t.get_text().replace("\n", " ").strip()
            if not t.get_visible() or not label:
                continue
            ext, poly = _ink(t, renderer)
            if ext.width < 0.5 or ext.height < 0.5:
                continue
            items.append((ax, label, ext, poly, t))
    for t in fig.texts:
        label = t.get_text().replace("\n", " ").strip()
        if not t.get_visible() or not label:
            continue
        ext, poly = _ink(t, renderer)
        items.append((None, label, ext, poly, t))
    return items


def _inside_anchor(text, ax):
    try:
        xy = text.get_transform().transform(text.get_unitless_position())
        frac = ax.transAxes.inverted().transform(xy)
    except Exception:
        return False
    return -0.02 <= frac[0] <= 1.02 and -0.02 <= frac[1] <= 1.02


def _data_boxes(ax, renderer):
    boxes = []
    ab = ax.get_window_extent(renderer)
    for coll in ax.collections:
        if not coll.get_visible():
            continue
        offs = coll.get_offsets()
        if offs is None or len(np.atleast_1d(offs)) == 0:
            continue
        xy = np.atleast_2d(coll.get_offset_transform().transform(np.atleast_2d(offs)))
        for x, y in xy:
            boxes.append(Bbox.from_bounds(x - 2.5, y - 2.5, 5, 5))
    for patch in ax.patches:
        if not patch.get_visible() or not ax.axison:
            continue
        bb = patch.get_window_extent(renderer)
        if bb.width < 1.5 or bb.height < 1.5:
            continue
        if bb.width > 0.92 * ab.width and bb.height > 0.92 * ab.height:
            continue
        boxes.append(bb)
    for line in ax.lines:
        if not line.get_visible():
            continue
        x = np.atleast_1d(line.get_xdata(orig=False))
        y = np.atleast_1d(line.get_ydata(orig=False))
        if len(x) < 2 or len(x) != len(y):
            continue
        xy = line.get_transform().transform(np.column_stack([np.asarray(x, float), np.asarray(y, float)]))
        span = np.ptp(xy, axis=0)
        if span[0] > 0.55 * ab.width or span[1] > 0.55 * ab.height:
            continue
        for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:]):
            boxes.append(Bbox.from_extents(min(x0, x1) - 1, min(y0, y1) - 1, max(x0, x1) + 1, max(y0, y1) + 1))
    return boxes


def _is_cell_number(name, text, ax):
    if name != "Fig6" or ax is None:
        return False
    label = text.get_text().strip()
    if not label.replace(".", "", 1).isdigit():
        return False
    return any(getattr(im, "get_array", lambda: None)() is not None for im in ax.images)


def audit(fig, name):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    items = _collect(fig, renderer)
    ticks = set()
    for ax in fig.axes:
        for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            ticks.add(id(t))
    fig_bb = fig.bbox
    rows = []
    for i, (ax, label, ext, poly, text) in enumerate(items):
        if ext.x0 < fig_bb.x0 - 3 or ext.y0 < fig_bb.y0 - 3 or ext.x1 > fig_bb.x1 + 3 or ext.y1 > fig_bb.y1 + 3:
            rows.append(dict(figure=name, kind="canvas", text=label,
                             overflow=(round(fig_bb.x0 - ext.x0, 1), round(fig_bb.y0 - ext.y0, 1),
                                       round(ext.x1 - fig_bb.x1, 1), round(ext.y1 - fig_bb.y1, 1))))
        for other in fig.axes:
            if other is ax or not other.axison:
                continue
            ob = other.get_window_extent(renderer)
            if ob.width > 0.9 * fig_bb.width and ob.height > 0.9 * fig_bb.height:
                continue
            area = _bb_area(ext, ob)
            if area >= 8:
                rows.append(dict(figure=name, kind="other_axes", text=label))
                break
        if ax is not None and id(text) not in ticks and text not in (ax.title, ax.xaxis.label, ax.yaxis.label):
            if _inside_anchor(text, ax):
                ab = ax.get_window_extent(renderer)
                if ext.x0 < ab.x0 - 1.5 or ext.x1 > ab.x1 + 1.5 or ext.y0 < ab.y0 - 1.5 or ext.y1 > ab.y1 + 1.5:
                    rows.append(dict(figure=name, kind="spine", text=label))
        if not _is_cell_number(name, text, ax) and ax is not None and id(text) not in ticks:
            for box in _data_boxes(ax, renderer):
                if _bb_area(ext, box) >= 6:
                    rows.append(dict(figure=name, kind="data", text=label))
                    break
        for j in range(i + 1, len(items)):
            ax2, label2, ext2, poly2, text2 = items[j]
            if text is text2:
                continue
            if id(text) in ticks and id(text2) in ticks:
                continue
            if _poly_overlap(poly, poly2):
                rows.append(dict(figure=name, kind="text", text=f"{label} | {label2}"))
    _tick_gaps(items, ticks, rows, name)
    for ax in fig.axes:
        legend = ax.get_legend()
        if legend is None:
            continue
        leg_bb = legend.get_window_extent(renderer)
        if leg_bb.x0 < fig_bb.x0 - 3 or leg_bb.y0 < fig_bb.y0 - 3 or leg_bb.x1 > fig_bb.x1 + 3 or leg_bb.y1 > fig_bb.y1 + 3:
            rows.append(dict(figure=name, kind="canvas", text="legend"))
        for other in fig.axes:
            if other is ax:
                continue
            if _bb_area(leg_bb, other.get_window_extent(renderer)) >= 8:
                rows.append(dict(figure=name, kind="other_axes", text="legend"))
                break
        for box in _data_boxes(ax, renderer):
            if _bb_area(leg_bb, box) >= 8:
                rows.append(dict(figure=name, kind="data", text="legend"))
                break
        for t in fig.texts:
            letter = t.get_text().strip()
            if not letter:
                continue
            if _bb_area(leg_bb, t.get_window_extent(renderer)) >= 4:
                rows.append(dict(figure=name, kind="legend_letter", text=f"legend | {letter}"))
    _line_hits(fig, items, ticks, rows, name, renderer)
    return rows


def _tick_gaps(items, ticks, rows, name):
    """Same-axis tick labels that touch, including a gap under 2 px."""
    groups = {}
    for ax, label, ext, poly, text in items:
        if ax is None or id(text) not in ticks or abs(text.get_rotation()) > 1:
            continue
        groups.setdefault(id(ax), []).append((label, ext))
    for group in groups.values():
        ordered = sorted(group, key=lambda r: (r[1].x0 + r[1].x1))
        for (lab_a, a), (lab_b, b) in zip(ordered, ordered[1:]):
            # side by side, not stacked
            if abs((a.y0 + a.y1) - (b.y0 + b.y1)) > 4:
                continue
            gap = b.x0 - a.x1
            if gap < 2.0:
                rows.append(dict(figure=name, kind="tick", text=f"{lab_a} | {lab_b}"))


def _line_hits(fig, items, ticks, rows, name, renderer):
    """Error-bar caps and vertical lines that pass through text."""
    from matplotlib.collections import LineCollection
    for ax in fig.axes:
        segments = []
        for line in ax.lines:
            if not line.get_visible():
                continue
            x = np.atleast_1d(np.asarray(line.get_xdata(orig=False), float))
            y = np.atleast_1d(np.asarray(line.get_ydata(orig=False), float))
            if len(x) < 2 or len(x) != len(y):
                continue
            xy = line.get_transform().transform(np.column_stack([x, y]))
            for p0, p1 in zip(xy[:-1], xy[1:]):
                segments.append((p0, p1))
        for coll in ax.collections:
            if not coll.get_visible() or not isinstance(coll, LineCollection):
                continue
            for seg in coll.get_segments():
                if len(seg) < 2:
                    continue
                xy = coll.get_transform().transform(np.asarray(seg, float))
                for p0, p1 in zip(xy[:-1], xy[1:]):
                    segments.append((p0, p1))
        for p0, p1 in segments:
            dx, dy = abs(p1[0] - p0[0]), abs(p1[1] - p0[1])
            vertical = dy > 8 and dy > dx * 3
            short = dx < 12 and dy < 12
            if not vertical and not short:
                continue
            kind = "vline" if vertical else "cap"
            for ax_t, label, ext, poly, text in items:
                if id(text) in ticks or ax_t is not ax:
                    continue
                if _segment_hits_poly(p0, p1, poly):
                    rows.append(dict(figure=name, kind=kind, text=label))
                    break
