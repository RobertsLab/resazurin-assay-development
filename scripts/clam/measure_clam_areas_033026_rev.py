#!/usr/bin/env python3
"""
Measure clam areas in 033026-images-rev plates using the red-outlined white
square as scale (1 cm × 1 cm in the object plane → 1 cm² reference area).
"""
from __future__ import annotations

import glob
import math
import os
import sys

import cv2
import numpy as np


def red_mask_hsv(bgr: np.ndarray) -> np.ndarray:
    hsv = cv2.cvtColor(bgr, cv2.COLOR_BGR2HSV)
    m1 = cv2.inRange(hsv, (0, 60, 60), (12, 255, 255))
    m2 = cv2.inRange(hsv, (168, 60, 60), (180, 255, 255))
    return cv2.bitwise_or(m1, m2)


def _select_calibration_red_contour(
    cnts: list,
    img_h: int,
    img_w: int,
) -> np.ndarray:
    """Pick ~square 1 cm frame; ignore large irregular red (e.g. labels)."""
    img_area = float(img_h * img_w)
    best: np.ndarray | None = None
    best_sc = -1e18
    for c in cnts:
        a = cv2.contourArea(c)
        if a < 350 or a > 0.015 * img_area:
            continue
        _, _, rw, rh = cv2.boundingRect(c)
        if rw < 28 or rh < 28:
            continue
        ar = rw / float(rh) if rh else 9.0
        if ar < 0.68 or ar > 1.48:
            continue
        side = 0.5 * (rw + rh)
        if side < 70 or side > 0.20 * min(img_h, img_w):
            continue
        sq = 1.0 - min(abs(1.0 - ar), 0.35)
        sz = 1.0 - min(abs(1.0 - side / 195.0), 0.6)
        sc = 0.55 * sq + 0.35 * sz + 0.10 * min(math.log(max(a, 10.0)) / 8.0, 1.0)
        if sc > best_sc:
            best_sc = sc
            best = c
    if best is None:
        return max(cnts, key=cv2.contourArea)
    return best


def calibration_scale_cm2_per_px(
    bgr: np.ndarray,
) -> tuple[float, dict, tuple[int, int, int, int]]:
    """
    Returns cm² per pixel (so area_cm2 = pixel_area * factor).
    Uses inner bright (“white”) region inside the red frame when possible.
    """
    h, w = bgr.shape[:2]
    red = red_mask_hsv(bgr)
    red = cv2.morphologyEx(red, cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8))
    red = cv2.morphologyEx(red, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    cnts, _ = cv2.findContours(red, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts:
        raise RuntimeError("No red contour found for calibration square.")

    c = _select_calibration_red_contour(cnts, h, w)
    x, y, rw, rh = cv2.boundingRect(c)
    pad = max(4, int(0.05 * max(rw, rh)))
    x0, y0 = max(0, x - pad), max(0, y - pad)
    x1, y1 = min(w, x + rw + pad), min(h, y + rh + pad)
    roi = bgr[y0:y1, x0:x1]
    gray = cv2.cvtColor(roi, cv2.COLOR_BGR2GRAY)
    rroi = red[y0:y1, x0:x1]

    # White interior: bright and not red (red dilated slightly to kill edge pixels)
    rdt = cv2.dilate(rroi, np.ones((5, 5), np.uint8))
    white = (gray >= 175) & (rdt == 0)
    white_u8 = white.astype(np.uint8) * 255
    white_u8 = cv2.morphologyEx(white_u8, cv2.MORPH_CLOSE, np.ones((7, 7), np.uint8))
    wc, _ = cv2.findContours(white_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not wc:
        # Fallback: use bounding box of red frame; assume square side = mean(w,h)
        side_px = float((rw + rh) / 2.0)
        ref_px2 = side_px**2
        meta = {"method": "red_bbox", "side_px": side_px, "ref_px2": ref_px2}
        return 1.0 / ref_px2, meta, (x, y, rw, rh)

    cw = max(wc, key=cv2.contourArea)
    ref_px2 = float(cv2.contourArea(cw))
    if ref_px2 < 100:
        side_px = float((rw + rh) / 2.0)
        ref_px2 = side_px**2
        meta = {"method": "red_bbox_tiny_white", "side_px": side_px, "ref_px2": ref_px2}
        return 1.0 / ref_px2, meta, (x, y, rw, rh)

    meta = {"method": "white_interior", "ref_px2": ref_px2, "red_bbox": (x, y, rw, rh)}
    return 1.0 / ref_px2, meta, (x, y, rw, rh)


def _exclude_rect_padded(
    exclude_xywh: tuple[int, int, int, int], pad: int, w: int, h: int
) -> tuple[int, int, int, int]:
    ex, ey, ew, eh = exclude_xywh
    X0, Y0 = max(0, ex - pad), max(0, ey - pad)
    X1, Y1 = min(w, ex + ew + pad), min(h, ey + eh + pad)
    return X0, Y0, X1, Y1


def _pick_central_clam_contour(
    binary: np.ndarray,
    cxl: float,
    cyl: float,
    rad: float,
    *,
    min_solidity: float = 0.34,
    max_centroid_dist_frac: float = 0.58,
) -> tuple[np.ndarray | None, float]:
    """
    Prefer a compact, roughly oval dark blob near the well center over edge rings.
    """
    cnts, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts:
        return None, 0.0

    def score_contour(c: np.ndarray, dist_cap: float) -> float | None:
        a = cv2.contourArea(c)
        if a < 80:
            return None
        M = cv2.moments(c)
        if M["m00"] < 1e-6:
            return None
        ccx = M["m10"] / M["m00"]
        ccy = M["m01"] / M["m00"]
        dist = math.hypot(ccx - cxl, ccy - cyl)
        if dist > dist_cap * rad:
            return None
        _, _, bw, bh = cv2.boundingRect(c)
        if bw < 8 or bh < 8:
            return None
        ar = bw / float(bh) if bh else 0.0
        if ar < 0.16 or ar > 6.2:
            return None
        hull = cv2.convexHull(c)
        hull_a = cv2.contourArea(hull)
        if hull_a > 1e-6 and (a / hull_a) < min_solidity:
            return None
        return a / (1.0 + 2.5 * (dist / max(rad, 1.0)) ** 2)

    best_c: np.ndarray | None = None
    best_score = -1.0
    for c in cnts:
        s = score_contour(c, max_centroid_dist_frac)
        if s is not None and s > best_score:
            best_score = s
            best_c = c
    if best_c is None:
        for c in cnts:
            s = score_contour(c, min(0.82, max_centroid_dist_frac + 0.22))
            if s is not None and s > best_score:
                best_score = s
                best_c = c
    if best_c is None:
        return None, 0.0
    return best_c, float(cv2.contourArea(best_c))


def _fallback_dark_blob(
    binary: np.ndarray,
    cxl: float,
    cyl: float,
    rad: float,
    amin: float,
    amax: float,
    *,
    dist_frac: float = 0.84,
) -> tuple[np.ndarray | None, float]:
    """If strict shape filters miss, take the largest plausible dark blob near center."""
    cnts, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    best: tuple[float, np.ndarray] | None = None
    for c in cnts:
        a = cv2.contourArea(c)
        if a < max(80.0, 0.5 * amin) or a > amax:
            continue
        M = cv2.moments(c)
        if M["m00"] < 1e-6:
            continue
        ccx = M["m10"] / M["m00"]
        ccy = M["m01"] / M["m00"]
        if math.hypot(ccx - cxl, ccy - cyl) > dist_frac * rad:
            continue
        if best is None or a > best[0]:
            best = (a, c)
    if best is None:
        return None, 0.0
    return best[1], float(cv2.contourArea(best[1]))


def _plate_bounds_from_luminance(
    blur: np.ndarray, H: int, W: int
) -> tuple[int, int, int, int]:
    """
    Infer plate crop from LAB L: styrofoam / margins stay very bright (~255);
    wells and clams are darker. Keeps fixed fractions as fallback.
    """
    y0s, y1s = int(0.10 * H), int(0.92 * H)
    x_slice = blur[y0s:y1s, :]
    col_med = np.median(x_slice, axis=0)
    dark_x = col_med < 251
    xs = np.flatnonzero(dark_x)
    if xs.size < max(100, int(0.04 * W)):
        return int(0.055 * W), int(0.945 * W), int(0.12 * H), int(0.91 * H)
    lx = max(0, int(xs.min()) - int(0.015 * W))
    rx = min(W, int(xs.max()) + int(0.015 * W) + 1)

    band = blur[y0s:y1s, lx:rx]
    row_med = np.median(band, axis=1)
    dark_y = row_med < 251
    ys = np.flatnonzero(dark_y)
    if ys.size < max(40, int(0.03 * H)):
        ty, by = int(0.12 * H), int(0.91 * H)
    else:
        ty = max(0, y0s + int(ys.min()) - int(0.012 * H))
        by = min(H, y0s + int(ys.max()) + int(0.012 * H) + 1)
    return lx, rx, ty, by


def measure_clams(
    bgr: np.ndarray,
    cm2_per_px: float,
    exclude_xywh: tuple[int, int, int, int],
    *,
    grid_rows: int = 6,
    grid_cols: int = 5,
) -> list[dict]:
    """
    One measurement per well: segment dark oval clams on a light well floor
    using Otsu on LAB L inside a circular well mask, then pick the central blob.
    """
    lab = cv2.cvtColor(bgr, cv2.COLOR_BGR2LAB)
    L = lab[:, :, 0]
    blur_full = cv2.GaussianBlur(L, (5, 5), 0)
    H, W = bgr.shape[:2]
    ex, ey, ew, eh = exclude_xywh
    cal_cx = ex + 0.5 * ew
    cal_cy = ey + 0.5 * eh

    lx, rx, ty, by = _plate_bounds_from_luminance(blur_full, H, W)
    cell_h = (by - ty) / float(grid_rows)
    cell_w = (rx - lx) / float(grid_cols)
    shrink = 0.12

    ref_px2 = 1.0 / cm2_per_px
    amin = 0.052 * ref_px2
    amax = 12.0 * ref_px2

    cells: list[tuple[int, int, int, int, int, int, float, float]] = []
    for gr in range(grid_rows):
        for gc in range(grid_cols):
            x0f = lx + (gc + shrink) * cell_w
            x1f = lx + (gc + 1 - shrink) * cell_w
            y0f = ty + (gr + shrink) * cell_h
            y1f = ty + (gr + 1 - shrink) * cell_h
            xi0, xi1 = int(x0f), int(x1f)
            yi0, yi1 = int(y0f), int(y1f)
            if xi1 - xi0 < 20 or yi1 - yi0 < 20:
                continue
            ccx = 0.5 * (xi0 + xi1)
            ccy = 0.5 * (yi0 + yi1)
            cells.append((gr, gc, xi0, xi1, yi0, yi1, ccx, ccy))

    skip_rc: tuple[int, int] | None = None
    for gr, gc, xi0, xi1, yi0, yi1, _, _ in cells:
        if xi0 <= cal_cx < xi1 and yi0 <= cal_cy < yi1:
            skip_rc = (gr, gc)
            break
    if skip_rc is None and cells:
        gr0, gc0, _, _, _, _, _, _ = min(
            cells,
            key=lambda t: (t[6] - cal_cx) ** 2 + (t[7] - cal_cy) ** 2,
        )
        skip_rc = (gr0, gc0)

    out: list[dict] = []
    for gr in range(grid_rows):
        for gc in range(grid_cols):
            if skip_rc is not None and (gr, gc) == skip_rc:
                continue
            x0f = lx + (gc + shrink) * cell_w
            x1f = lx + (gc + 1 - shrink) * cell_w
            y0f = ty + (gr + shrink) * cell_h
            y1f = ty + (gr + 1 - shrink) * cell_h
            xi0, xi1 = int(x0f), int(x1f)
            yi0, yi1 = int(y0f), int(y1f)
            if xi1 - xi0 < 20 or yi1 - yi0 < 20:
                continue

            roi = blur_full[yi0:yi1, xi0:xi1]
            if roi.size == 0:
                continue
            rh, rw = roi.shape
            cyl, cxl = rh * 0.5, rw * 0.5
            yy, xx = np.ogrid[:rh, :rw]
            rad = 0.44 * float(min(rh, rw))
            well_disk = (yy - cyl) ** 2 + (xx - cxl) ** 2 <= rad**2
            samp = roi[well_disk]
            if samp.size < 80:
                continue

            p_lo = float(np.percentile(samp, 18))
            p_mid = float(np.percentile(samp, 45))
            p_hi = float(np.percentile(samp, 80))
            spread = p_hi - p_lo
            if spread < 3.0 or (spread < 4.0 and not (gc == 0 or gc == grid_cols - 1)):
                continue

            col = samp.reshape(-1, 1).astype(np.uint8)
            t_otsu_i, _ = cv2.threshold(col, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
            t_otsu = float(t_otsu_i)

            t_pct = 0.38 * p_hi + 0.62 * p_lo
            T_base = min(t_otsu, t_pct + 4.0)
            edge_slack = min(11, max(3, 0.04 * (p_hi - p_lo)))
            cuts = [
                T_base + edge_slack,
                T_base + edge_slack + 12,
                T_base + edge_slack + 22,
                max(0.0, T_base + edge_slack - 10),
                0.5 * p_mid + 0.5 * p_lo + 3.0,
                0.44 * p_hi + 0.56 * p_lo,
                0.48 * p_hi + 0.52 * p_lo + 2.0,
                0.5 * (p_lo + p_hi) + 4.0,
                t_otsu + edge_slack + 6.0,
            ]

            edge_col = gc == 0 or gc == grid_cols - 1
            dist_frac = 0.82 if edge_col else 0.70

            c = None
            a = 0.0
            for t_cut in cuts:
                dark = (roi.astype(np.float32) < float(t_cut)) & well_disk
                binary = dark.astype(np.uint8) * 255
                if int(binary.sum()) < 100:
                    continue
                binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
                binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, np.ones((11, 11), np.uint8))
                c_try, a_try = _pick_central_clam_contour(
                    binary,
                    cxl,
                    cyl,
                    rad,
                    min_solidity=0.24,
                    max_centroid_dist_frac=dist_frac,
                )
                if c_try is None or a_try < 1.0:
                    c_try, a_try = _fallback_dark_blob(
                        binary, cxl, cyl, rad, amin, amax, dist_frac=dist_frac + 0.08
                    )
                if c_try is not None and amin <= a_try <= amax:
                    c, a = c_try, a_try
                    break
            if c is None or a < 1.0:
                continue
            M = cv2.moments(c)
            if M["m00"] == 0:
                continue
            cx = M["m10"] / M["m00"] + xi0
            cy = M["m01"] / M["m00"] + yi0
            x, y, brw, brh = cv2.boundingRect(c)
            x += xi0
            y += yi0
            cnt_g = c.astype(np.int32).copy()
            cnt_g[:, 0, 0] += xi0
            cnt_g[:, 0, 1] += yi0
            out.append(
                {
                    "cx": cx,
                    "cy": cy,
                    "area_px": float(a),
                    "area_cm2": float(a * cm2_per_px),
                    "bbox": (x, y, brw, brh),
                    "grid_r": gr + 1,
                    "grid_c": gc + 1,
                    "contour": cnt_g,
                }
            )
    return out


def format_report(
    image_name: str,
    cm2_per_px: float,
    meta: dict,
    clams: list[dict],
) -> list[str]:
    lines: list[str] = []
    lines.append("=" * 72)
    lines.append(f"  Image: {image_name}")
    lines.append(
        f"  Scale: 1 cm² reference = {1.0 / cm2_per_px:.1f} px²  |  calibration: {meta}"
    )
    lines.append(f"  Clams detected: {len(clams)}")
    lines.append("")
    lines.append("  Idx  RowCol  Center(x,y)px   Area_px   Area_cm²")
    lines.append("  ---  ------  -------------   -------   --------")
    for i, q in enumerate(clams, start=1):
        cx, cy = q["cx"], q["cy"]
        pos = q.get("pos", "  ?  ")
        lines.append(
            f"  {i:3d}  {pos:6s}  ({cx:7.1f},{cy:7.1f})   {q['area_px']:7.1f}   {q['area_cm2']:8.3f}"
        )
    if clams:
        areas = [q["area_cm2"] for q in clams]
        mean = sum(areas) / len(areas)
        var = sum((a - mean) ** 2 for a in areas) / len(areas)
        sd = var**0.5
        lines.append("")
        lines.append(
            f"  Summary: n={len(clams)}  mean={mean:.3f} cm²  SD={sd:.3f}  "
            f"min={min(areas):.3f}  max={max(areas):.3f}"
        )
    lines.append("")
    return lines


def draw_clam_overlays(bgr: np.ndarray, clams: list[dict]) -> np.ndarray:
    """Draw segmentation contours and 1-based idx (same order as clam_area_estimates.txt)."""
    vis = bgr.copy()
    h = bgr.shape[0]
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_scale = max(0.7, min(2.2, h / 1400.0))
    thickness = max(2, int(round(font_scale * 1.5)))
    outline = thickness + 2
    green = (60, 200, 80)

    for i, q in enumerate(clams, start=1):
        cnt = q.get("contour")
        if cnt is not None and len(cnt) >= 3:
            cv2.drawContours(vis, [cnt], -1, green, 2, lineType=cv2.LINE_AA)

        cx, cy = int(round(q["cx"])), int(round(q["cy"]))
        label = str(i)
        (tw, th), baseline = cv2.getTextSize(label, font, font_scale, thickness)
        tx = max(0, min(vis.shape[1] - tw, cx - tw // 2))
        ty = max(th + 4, cy - 12)

        cv2.putText(
            vis,
            label,
            (tx, ty),
            font,
            font_scale,
            (0, 0, 0),
            outline,
            lineType=cv2.LINE_AA,
        )
        cv2.putText(
            vis,
            label,
            (tx, ty),
            font,
            font_scale,
            (255, 255, 255),
            thickness,
            lineType=cv2.LINE_AA,
        )
    return vis


def assign_row_col(clams: list[dict], n_rows: int = 6, n_cols: int = 5) -> None:
    """Assign R{C} labels (prefer plate grid index from segmentation)."""
    if not clams:
        return
    infer: list[dict] = []
    for q in clams:
        if "grid_r" in q and "grid_c" in q:
            q["row"], q["col"] = q["grid_r"], q["grid_c"]
            q["pos"] = f"R{q['row']}C{q['col']}"
        else:
            infer.append(q)
    if not infer:
        return
    infer.sort(key=lambda q: (q["cy"], q["cx"]))
    ys = [q["cy"] for q in infer]
    ymin, ymax = min(ys), max(ys)
    xs = [q["cx"] for q in infer]
    xmin, xmax = min(xs), max(xs)

    def row_of(cy: float) -> int:
        if ymax == ymin:
            return 1
        t = (cy - ymin) / (ymax - ymin)
        r = int(round(1 + t * (n_rows - 1)))
        return max(1, min(n_rows, r))

    def col_of(cx: float) -> int:
        if xmax == xmin:
            return 1
        t = (cx - xmin) / (xmax - xmin)
        c = int(round(1 + t * (n_cols - 1)))
        return max(1, min(n_cols, c))

    for q in infer:
        q["row"] = row_of(q["cy"])
        q["col"] = col_of(q["cx"])
        q["pos"] = f"R{q['row']}C{q['col']}"


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.abspath(os.path.join(here, "..", ".."))
    img_dir = os.path.join(repo, "data", "clam", "033026-images-rev")
    out_path = os.path.join(img_dir, "clam_area_estimates.txt")
    paths = sorted(
        p
        for p in glob.glob(os.path.join(img_dir, "*.png"))
        if "annotated" not in os.path.basename(p).lower()
    )
    if not paths:
        print(f"No PNG files in {img_dir}", file=sys.stderr)
        return 1

    header = [
        "=" * 72,
        "  CLAM AREAS — 033026-images-rev",
        "  Reference: white square, red outline = 1 cm cube face in the image plane (1 cm²).",
        "  Method: auto plate crop from L; ~square red frame pick; 6×5 grid; per-well "
        "LAB L + Otsu/percentile cuts; circular ROI; central compact blob; area → cm².",
        "=" * 72,
        "",
    ]
    all_lines = header

    for p in paths:
        bgr = cv2.imread(p)
        if bgr is None:
            print(f"Failed to read {p}", file=sys.stderr)
            continue
        try:
            cm2_per_px, meta, excl = calibration_scale_cm2_per_px(bgr)
        except RuntimeError as e:
            print(f"{os.path.basename(p)}: {e}", file=sys.stderr)
            continue
        clams = measure_clams(bgr, cm2_per_px, excl)
        assign_row_col(clams)
        clams.sort(key=lambda q: (q["row"], q["col"], q["cy"], q["cx"]))
        name = os.path.basename(p)
        all_lines.extend(format_report(name, cm2_per_px, meta, clams))

        ann_dir = os.path.join(img_dir, "annotated")
        os.makedirs(ann_dir, exist_ok=True)
        stem, _ = os.path.splitext(name)
        ann_path = os.path.join(ann_dir, f"{stem}_annotated.png")
        overlay = draw_clam_overlays(bgr, clams)
        if not cv2.imwrite(ann_path, overlay):
            print(f"Failed to write {ann_path}", file=sys.stderr)
        else:
            print(f"Wrote {ann_path}")

    all_lines.append("=" * 72)
    all_lines.append(
        "  Notes: 6×5 well layout (30 positions). The well covered by the scale square "
        "is omitted. Reported area is the segmented dark region in each remaining well "
        "(planar cm² from the photo), not a 3D shell volume."
    )
    all_lines.append(
        "  QC images: segmented contour + idx (same numbering as this file) in "
        "subfolder annotated/*_annotated.png"
    )
    all_lines.append("=" * 72)

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(all_lines))
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
