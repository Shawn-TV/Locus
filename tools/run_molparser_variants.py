#!/usr/bin/env python3
"""Run MolParser on safer preprocessed image variants for weak recognitions."""

from __future__ import annotations

import argparse
import base64
import csv
import io
import json
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageChops, ImageFilter, ImageOps

try:
    import cv2
    import numpy as np
except Exception:  # pragma: no cover - color cleanup is optional.
    cv2 = None
    np = None


API_URL = "https://ocsr.dp.tech/mol/img2mol"
FIELDS = [
    "file_name",
    "caption",
    "smi",
    "score",
    "markush",
    "n_results",
    "bbox",
    "trace_id",
    "error",
    "seconds",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def successful_rows(path: Path) -> dict[str, tuple[str, float]]:
    best: dict[str, tuple[str, float]] = {}
    for row in read_csv(path):
        if row.get("error") or not row.get("file_name"):
            continue
        caption = row.get("caption") or (f"{row.get('smi', '')}<sep>" if row.get("smi") else "")
        if not caption:
            continue
        try:
            score = float(row.get("score") or 0.0)
        except ValueError:
            score = 0.0
        old = best.get(row["file_name"])
        if old is None or score > old[1]:
            best[row["file_name"]] = (caption, score)
    return best


def names_from_args(args: argparse.Namespace) -> list[str]:
    names: list[str] = []
    if args.names:
        names.extend(args.names)
    if args.names_file:
        names.extend(
            line.strip()
            for line in args.names_file.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    if args.from_csv:
        best = successful_rows(args.from_csv)
        for path in sorted(args.input_dir.iterdir()):
            if path.suffix.lower() not in {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff"}:
                continue
            if path.name not in best or best[path.name][1] < args.threshold:
                names.append(path.name)
    return sorted(dict.fromkeys(names))


def content_bbox(img: Image.Image, threshold: int) -> tuple[int, int, int, int] | None:
    rgb = img.convert("RGB")
    bg = Image.new("RGB", rgb.size, "white")
    diff = ImageChops.difference(rgb, bg).convert("L")
    return diff.point(lambda px: 255 if px > threshold else 0).getbbox()


def padded(img: Image.Image, pad: int) -> Image.Image:
    canvas = Image.new("RGB", (img.width + 2 * pad, img.height + 2 * pad), "white")
    canvas.paste(img.convert("RGB"), (pad, pad))
    return canvas


def threshold_ink(img: Image.Image, cutoff: int = 210) -> Image.Image:
    gray = ImageOps.grayscale(img.convert("RGB"))
    bw = gray.point(lambda px: 0 if px < cutoff else 255, mode="1").convert("RGB")
    return bw.filter(ImageFilter.MedianFilter(size=3))


def color_foreground(img: Image.Image) -> Image.Image | None:
    if cv2 is None or np is None:
        return None
    rgb = np.array(img.convert("RGB"))
    maxc = rgb.max(axis=2).astype("int16")
    minc = rgb.min(axis=2).astype("int16")
    luma = (0.299 * rgb[:, :, 0] + 0.587 * rgb[:, :, 1] + 0.114 * rgb[:, :, 2]).astype("uint8")
    saturation = maxc - minc

    colored = ((saturation > 45) & (luma < 245)).astype("uint8")
    if int(colored.sum()) < 80:
        return None
    dark = (luma < 105).astype("uint8")
    kernel = np.ones((25, 25), dtype="uint8")
    near_colored = cv2.dilate(colored, kernel, iterations=1)
    mask = ((colored > 0) | ((dark > 0) & (near_colored > 0))).astype("uint8")

    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask, connectivity=8)
    keep = np.zeros_like(mask)
    for label in range(1, num_labels):
        area = stats[label, cv2.CC_STAT_AREA]
        component = labels == label
        if area >= 18 and (colored[component].sum() > 0 or area >= 120):
            keep[component] = 1
    if int(keep.sum()) < 80:
        return None

    out = np.full(rgb.shape, 255, dtype="uint8")
    out[keep > 0] = 0
    ys, xs = np.where(keep > 0)
    x0, x1 = max(0, xs.min() - 8), min(out.shape[1], xs.max() + 9)
    y0, y1 = max(0, ys.min() - 8), min(out.shape[0], ys.max() + 9)
    return Image.fromarray(out[y0:y1, x0:x1])


def variants(path: Path, include_aggressive: bool) -> Iterable[tuple[str, Image.Image]]:
    img = Image.open(path).convert("RGB")
    yield "border80", padded(img, 80)

    bbox = content_bbox(img, threshold=12)
    cropped = img.crop(bbox) if bbox else img
    for pad in (24, 60, 120):
        yield f"crop_pad{pad}", padded(cropped, pad)

    up2 = cropped.resize((cropped.width * 2, cropped.height * 2), Image.Resampling.LANCZOS)
    yield "crop_up2_pad60", padded(up2, 60)

    if include_aggressive:
        color_cleaned = color_foreground(img)
        if color_cleaned is not None:
            yield "color_fg_pad60", padded(color_cleaned, 60)
            up = color_cleaned.resize(
                (color_cleaned.width * 2, color_cleaned.height * 2),
                Image.Resampling.LANCZOS,
            )
            yield "color_fg_up2_pad60", padded(up, 60)
        for source_name, source in (("orig", img), ("crop", cropped)):
            bw = threshold_ink(source)
            yield f"bw_{source_name}_pad60", padded(bw, 60)
            up = bw.resize((bw.width * 2, bw.height * 2), Image.Resampling.LANCZOS)
            yield f"bw_{source_name}_up2_pad60", padded(up, 60)


def post_image(name: str, variant_name: str, img: Image.Image, timeout: float, retries: int) -> dict[str, object]:
    started = time.time()
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    payload = json.dumps({"base64_img": base64.b64encode(buf.getvalue()).decode("ascii")}).encode()
    headers = {
        "Content-Type": "application/json",
        "User-Agent": "Mozilla/5.0 Bohrium-V2-OCSR",
    }
    last_error = ""
    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(API_URL, data=payload, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                parsed = json.loads(response.read().decode("utf-8"))
            if parsed.get("code") != 0:
                raise RuntimeError(parsed.get("msg") or parsed.get("code"))
            results = parsed.get("data") or []
            best = results[0] if results else {}
            return {
                "file_name": name,
                "caption": best.get("caption", ""),
                "smi": best.get("smi", ""),
                "score": best.get("score", ""),
                "markush": best.get("markush", ""),
                "n_results": len(results),
                "bbox": json.dumps(best.get("bbox", ""), ensure_ascii=False),
                "trace_id": f"variant:{variant_name}:{parsed.get('trace_id', '')}",
                "error": "",
                "seconds": round(time.time() - started, 3),
            }
        except (urllib.error.URLError, TimeoutError, RuntimeError, json.JSONDecodeError) as exc:
            last_error = repr(exc)
            if attempt < retries:
                time.sleep(1.5 * (attempt + 1))
    return {
        "file_name": name,
        "caption": "",
        "smi": "",
        "score": "",
        "markush": "",
        "n_results": 0,
        "bbox": "",
        "trace_id": f"variant:{variant_name}",
        "error": last_error,
        "seconds": round(time.time() - started, 3),
    }


def append_rows(output: Path, rows: list[dict[str, object]]) -> None:
    exists = output.exists() and output.stat().st_size > 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if not exists:
            writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--from-csv", type=Path)
    parser.add_argument("--threshold", type=float, default=0.95)
    parser.add_argument("--names-file", type=Path)
    parser.add_argument("--names", nargs="*")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--retries", type=int, default=1)
    parser.add_argument("--include-aggressive", action="store_true")
    args = parser.parse_args()

    names = names_from_args(args)
    jobs: list[tuple[str, str, Image.Image]] = []
    for name in names:
        path = args.input_dir / name
        if path.exists():
            jobs.extend((name, variant_name, img) for variant_name, img in variants(path, args.include_aggressive))

    print(
        json.dumps(
            {"event": "start", "names": len(names), "jobs": len(jobs), "workers": args.workers},
            ensure_ascii=False,
        ),
        flush=True,
    )

    completed = 0
    success = 0
    rows: list[dict[str, object]] = []
    started = time.time()
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {
            pool.submit(post_image, name, variant_name, img, args.timeout, args.retries): (name, variant_name)
            for name, variant_name, img in jobs
        }
        for future in as_completed(futures):
            row = future.result()
            completed += 1
            if not row.get("error") and (row.get("caption") or row.get("smi")):
                success += 1
                rows.append(row)
            if len(rows) >= 50:
                append_rows(args.output, rows)
                rows.clear()
            if completed % 50 == 0 or completed == len(jobs):
                print(
                    json.dumps(
                        {
                            "event": "progress",
                            "completed": completed,
                            "jobs": len(jobs),
                            "success": success,
                            "elapsed": round(time.time() - started, 3),
                        },
                        ensure_ascii=False,
                    ),
                    flush=True,
                )
    if rows:
        append_rows(args.output, rows)
    print(
        json.dumps(
            {"event": "complete", "completed": completed, "success": success, "elapsed": round(time.time() - started, 3)},
            ensure_ascii=False,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
