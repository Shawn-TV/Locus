#!/usr/bin/env python3
"""Audit V1 predictions by rendering E-SMILES and ranking risky rows.

This script is intentionally conservative.  Render/image similarity is only a
triage signal because a correct molecule can be drawn in many layouts.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import re
import textwrap
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from skimage.metrics import structural_similarity


RDLogger.DisableLog("rdApp.*")


ANNOT_RE = re.compile(r"<([ar])>(\d+):([^<]+)</[ar]>")


@dataclass
class Candidate:
    source: str
    caption: str
    score: float
    valid: bool
    canonical: str
    atoms: int


def read_csv(path: Path, encoding: str = "utf-8") -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding=encoding, newline="") as f:
        return list(csv.DictReader(f))


def ensure_esmi(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return ""
    return value if "<sep>" in value else f"{value}<sep>"


def base_part(esmi: str) -> str:
    return ensure_esmi(esmi).split("<sep>")[0]


def annotations(esmi: str) -> dict[int, str]:
    labels: dict[int, str] = {}
    for _kind, idx, label in ANNOT_RE.findall(ensure_esmi(esmi)):
        labels[int(idx)] = label
    return labels


def atom_annotations(esmi: str) -> dict[int, str]:
    labels: dict[int, str] = {}
    for kind, idx, label in ANNOT_RE.findall(ensure_esmi(esmi)):
        if kind == "a":
            labels[int(idx)] = label
    return labels


def ring_annotations(esmi: str) -> dict[int, str]:
    labels: dict[int, str] = {}
    for kind, idx, label in ANNOT_RE.findall(ensure_esmi(esmi)):
        if kind == "r":
            labels[int(idx)] = label
    return labels


def mol_from_esmi(esmi: str) -> Chem.Mol | None:
    base = base_part(esmi)
    if not base:
        return None
    try:
        return Chem.MolFromSmiles(base)
    except Exception:
        return None


def canonical_base(esmi: str) -> tuple[bool, str, int]:
    mol = mol_from_esmi(esmi)
    if mol is None:
        return False, "", 0
    return True, Chem.MolToSmiles(mol, isomericSmiles=True), mol.GetNumAtoms()


def render_esmi(esmi: str, out_path: Path, size: tuple[int, int] = (420, 320)) -> bool:
    mol = mol_from_esmi(esmi)
    if mol is None:
        return False
    mol = Chem.Mol(mol)
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    options = drawer.drawOptions()
    options.clearBackground = True
    options.bondLineWidth = 2.0
    for idx, label in atom_annotations(esmi).items():
        if 0 <= idx < mol.GetNumAtoms():
            options.atomLabels[idx] = label
    drawer.DrawMolecule(mol)
    ring_label_positions: list[tuple[float, float, str]] = []
    for ring_idx, label in ring_annotations(esmi).items():
        rings = mol.GetRingInfo().AtomRings()
        if 0 <= ring_idx < len(rings) and rings[ring_idx]:
            points = [drawer.GetDrawCoords(atom_idx) for atom_idx in rings[ring_idx]]
            x = sum(point.x for point in points) / len(points)
            y = sum(point.y for point in points) / len(points)
            ring_label_positions.append((x - 14, y - 8, label))
    drawer.FinishDrawing()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    drawing = drawer.GetDrawingText()
    if ring_label_positions:
        img = Image.open(io.BytesIO(drawing)).convert("RGB")
        draw = ImageDraw.Draw(img)
        try:
            font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 20)
        except Exception:
            font = None
        for x, y, label in ring_label_positions:
            draw.text((x, y), label, fill=(180, 0, 0), font=font)
        img.save(out_path)
    else:
        out_path.write_bytes(drawing)
    return True


def content_bbox(img: Image.Image, threshold: int = 18) -> tuple[int, int, int, int] | None:
    rgb = img.convert("RGB")
    bg = Image.new("RGB", rgb.size, "white")
    diff = ImageChops.difference(rgb, bg).convert("L")
    return diff.point(lambda px: 255 if px > threshold else 0).getbbox()


def normalize_for_compare(path: Path | None, img: Image.Image | None = None) -> np.ndarray | None:
    if img is None:
        if path is None or not path.exists():
            return None
        img = Image.open(path).convert("RGB")
    bbox = content_bbox(img, 14)
    if bbox:
        img = img.crop(bbox)
    img = ImageOps.pad(img, (256, 256), color="white", method=Image.Resampling.LANCZOS)
    gray = ImageOps.grayscale(img)
    arr = np.array(gray)
    arr = np.where(arr < 230, 0, 255).astype("uint8")
    return arr


def phash(arr: np.ndarray) -> np.ndarray:
    small = cv2.resize(arr, (32, 32), interpolation=cv2.INTER_AREA).astype("float32")
    dct = cv2.dct(small)
    block = dct[:8, :8].copy()
    vals = block.flatten()
    median = np.median(vals[1:])
    return vals > median


def compare_images(original_path: Path, rendered_path: Path) -> tuple[float, float, float]:
    a = normalize_for_compare(original_path)
    b = normalize_for_compare(rendered_path)
    if a is None or b is None:
        return math.nan, math.nan, math.nan
    hamming = float(np.count_nonzero(phash(a) != phash(b)))
    try:
        ssim = float(structural_similarity(a, b, data_range=255))
    except Exception:
        ssim = math.nan
    ink_a = float((a == 0).mean())
    ink_b = float((b == 0).mean())
    return hamming, ssim, abs(ink_a - ink_b)


def image_features(path: Path) -> tuple[float, float, int, int]:
    img = Image.open(path).convert("RGB")
    arr = np.array(img).astype("int16")
    saturation = arr.max(axis=2) - arr.min(axis=2)
    luma = (0.299 * arr[:, :, 0] + 0.587 * arr[:, :, 1] + 0.114 * arr[:, :, 2])
    colored = ((saturation > 45) & (luma < 245)).mean()
    ink = (luma < 230).mean()
    return float(colored), float(ink), img.width, img.height


def load_v1(zip_path: Path) -> dict[str, str]:
    with zipfile.ZipFile(zip_path) as zf:
        rows = csv.DictReader(zf.read("submission.csv").decode("utf-8-sig").splitlines())
        return {row["file_name"]: ensure_esmi(row["e_smiles"]) for row in rows}


def load_molparser(path: Path, source: str) -> dict[str, list[Candidate]]:
    out: dict[str, list[Candidate]] = defaultdict(list)
    for row in read_csv(path):
        name = row.get("file_name", "")
        if not name or row.get("error"):
            continue
        caption = ensure_esmi(row.get("caption", "")) or ensure_esmi(row.get("smi", ""))
        if not caption:
            continue
        try:
            score = float(row.get("score") or 0.0)
        except ValueError:
            score = 0.0
        valid, canon, atoms = canonical_base(caption)
        out[name].append(Candidate(source, caption, score, valid, canon, atoms))
    return out


def load_molscribe(path: Path) -> dict[str, list[Candidate]]:
    out: dict[str, list[Candidate]] = defaultdict(list)
    for row in read_csv(path):
        name = row.get("file_name", "")
        if not name or row.get("error"):
            continue
        for key, source in (("e_smiles", "molscribe"), ("raw_smiles", "molscribe_raw")):
            caption = ensure_esmi(row.get(key, ""))
            if not caption or "<invalid>" in caption:
                continue
            valid, canon, atoms = canonical_base(caption)
            out[name].append(Candidate(source, caption, 0.0, valid, canon, atoms))
    return out


def merge_candidate_maps(*maps: dict[str, list[Candidate]]) -> dict[str, list[Candidate]]:
    merged: dict[str, list[Candidate]] = defaultdict(list)
    for mapping in maps:
        for name, candidates in mapping.items():
            merged[name].extend(candidates)
    return merged


def best_candidate(candidates: list[Candidate]) -> Candidate | None:
    if not candidates:
        return None
    return max(candidates, key=lambda c: (c.valid, c.score, c.atoms, len(c.caption)))


def best_different_candidate(candidates: list[Candidate], v1_canon: str) -> Candidate | None:
    different = [c for c in candidates if c.valid and c.canonical and c.canonical != v1_canon]
    if not different:
        return None
    return max(different, key=lambda c: (c.score >= 0.95, c.score, c.atoms, len(c.caption)))


def score_row(
    name: str,
    v1: str,
    valid: bool,
    v1_canon: str,
    atoms: int,
    orig_candidates: list[Candidate],
    all_candidates: list[Candidate],
    phash_distance: float,
    ssim: float,
    ink_delta: float,
    colored_ratio: float,
) -> tuple[float, list[str], Candidate | None, Candidate | None]:
    reasons: list[str] = []
    risk = 0.0
    if not valid:
        risk += 120
        reasons.append("V1 base SMILES invalid")

    orig_best = best_candidate(orig_candidates)
    all_best = best_candidate(all_candidates)
    diff_best = best_different_candidate(all_candidates, v1_canon)

    if orig_best is None:
        risk += 36
        reasons.append("no original MolParser success")
    else:
        score = orig_best.score
        if score < 0.5:
            risk += 45
            reasons.append(f"low original score {score:.3f}")
        elif score < 0.7:
            risk += 34
            reasons.append(f"medium-low original score {score:.3f}")
        elif score < 0.85:
            risk += 24
            reasons.append(f"weak original score {score:.3f}")
        elif score < 0.95:
            risk += 14
            reasons.append(f"sub-0.95 original score {score:.3f}")

    if diff_best and diff_best.score >= 0.95:
        risk += 22
        reasons.append(f"high-score different candidate {diff_best.source}:{diff_best.score:.3f}")
    elif diff_best and diff_best.score >= 0.7:
        risk += 12
        reasons.append(f"different candidate {diff_best.source}:{diff_best.score:.3f}")

    if annotations(v1):
        risk += 5
        reasons.append("has E-SMILES labels")

    if colored_ratio > 0.01:
        risk += min(12, colored_ratio * 300)
        reasons.append(f"colored/label-heavy image {colored_ratio:.3f}")

    if not math.isnan(phash_distance) and phash_distance > 28:
        risk += min(18, (phash_distance - 28) * 0.8)
        reasons.append(f"render visual mismatch phash {phash_distance:.0f}")
    if not math.isnan(ssim) and ssim < 0.18:
        risk += 8
        reasons.append(f"low render SSIM {ssim:.2f}")
    if not math.isnan(ink_delta) and ink_delta > 0.08:
        risk += 5
        reasons.append(f"ink density delta {ink_delta:.2f}")

    if atoms <= 4:
        risk += 5
        reasons.append("very small/label-like structure")

    return risk, reasons, orig_best, all_best or orig_best


def short(value: str, limit: int = 100) -> str:
    value = value.replace("\n", " ")
    return value if len(value) <= limit else value[: limit - 1] + "…"


def make_review_sheets(
    rows: list[dict[str, str]],
    pic_dir: Path,
    render_dir: Path,
    out_dir: Path,
    top_n: int,
    per_page: int = 12,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    selected = rows[:top_n]
    tile_w, tile_h = 900, 300
    for page_start in range(0, len(selected), per_page):
        page = selected[page_start : page_start + per_page]
        sheet = Image.new("RGB", (tile_w, tile_h * len(page)), "white")
        for i, row in enumerate(page):
            tile = Image.new("RGB", (tile_w, tile_h), "white")
            draw = ImageDraw.Draw(tile)
            name = row["file_name"]
            orig = Image.open(pic_dir / name).convert("RGB")
            rendered_path = render_dir / name
            rendered = Image.open(rendered_path).convert("RGB") if rendered_path.exists() else Image.new("RGB", (420, 320), "white")
            for img, x in ((orig, 10), (rendered, 350)):
                img = ImageOps.contain(img, (320, 210), Image.Resampling.LANCZOS)
                frame = Image.new("RGB", (330, 220), "white")
                frame.paste(img, ((330 - img.width) // 2, (220 - img.height) // 2))
                tile.paste(frame, (x, 45))
                draw.rectangle((x, 45, x + 329, 264), outline=(210, 210, 210), width=1)
            header = f"{name} risk={float(row['risk']):.1f} score={row['orig_score']} source={row['best_source']}"
            draw.text((10, 8), header, fill="black")
            draw.text((10, 270), "orig", fill=(80, 80, 80))
            draw.text((350, 270), "V1 render", fill=(80, 80, 80))
            reason = "\n".join(textwrap.wrap(row["reasons"], 38))
            draw.multiline_text((690, 45), reason[:260], fill=(160, 0, 0), spacing=3)
            candidate = "best diff: " + short(row["diff_candidate"], 120)
            draw.multiline_text((690, 170), "\n".join(textwrap.wrap(candidate, 38))[:170], fill=(0, 70, 120), spacing=3)
            sheet.paste(tile, (0, i * tile_h))
        sheet.save(out_dir / f"review_page_{page_start // per_page + 1:02d}.png")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v1-zip", required=True, type=Path)
    parser.add_argument("--pic-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--molparser-csv", nargs="+", required=True, type=Path)
    parser.add_argument("--molscribe-csv", required=True, type=Path)
    parser.add_argument("--review-top", type=int, default=180)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    render_dir = args.out_dir / "rendered"
    v1 = load_v1(args.v1_zip)

    orig_map = load_molparser(args.molparser_csv[0], "molparser_original")
    other_maps = [load_molparser(path, f"molparser_{path.stem}") for path in args.molparser_csv[1:]]
    molscribe_map = load_molscribe(args.molscribe_csv)
    all_map = merge_candidate_maps(orig_map, *other_maps, molscribe_map)

    out_rows: list[dict[str, str]] = []
    for idx, (name, esmi) in enumerate(sorted(v1.items()), start=1):
        rendered_path = render_dir / name
        rendered_ok = render_esmi(esmi, rendered_path)
        valid, canon, atoms = canonical_base(esmi)
        if rendered_ok:
            phash_distance, ssim, ink_delta = compare_images(args.pic_dir / name, rendered_path)
        else:
            phash_distance, ssim, ink_delta = math.nan, math.nan, math.nan
        colored_ratio, ink_ratio, width, height = image_features(args.pic_dir / name)
        risk, reasons, orig_best, all_best = score_row(
            name,
            esmi,
            valid,
            canon,
            atoms,
            orig_map.get(name, []),
            all_map.get(name, []),
            phash_distance,
            ssim,
            ink_delta,
            colored_ratio,
        )
        diff = best_different_candidate(all_map.get(name, []), canon)
        out_rows.append(
            {
                "file_name": name,
                "risk": f"{risk:.3f}",
                "reasons": "; ".join(reasons),
                "v1_valid": str(valid),
                "v1_atoms": str(atoms),
                "orig_score": "" if orig_best is None else f"{orig_best.score:.6g}",
                "best_source": "" if all_best is None else all_best.source,
                "best_score": "" if all_best is None else f"{all_best.score:.6g}",
                "diff_source": "" if diff is None else diff.source,
                "diff_score": "" if diff is None else f"{diff.score:.6g}",
                "phash": "" if math.isnan(phash_distance) else f"{phash_distance:.3f}",
                "ssim": "" if math.isnan(ssim) else f"{ssim:.5f}",
                "ink_delta": "" if math.isnan(ink_delta) else f"{ink_delta:.5f}",
                "colored_ratio": f"{colored_ratio:.6f}",
                "ink_ratio": f"{ink_ratio:.6f}",
                "image_width": str(width),
                "image_height": str(height),
                "v1_esmi": esmi,
                "orig_candidate": "" if orig_best is None else orig_best.caption,
                "best_candidate": "" if all_best is None else all_best.caption,
                "diff_candidate": "" if diff is None else diff.caption,
            }
        )
        if idx % 250 == 0:
            print({"event": "progress", "rows": idx}, flush=True)

    out_rows.sort(key=lambda row: float(row["risk"]), reverse=True)
    audit_path = args.out_dir / "v1_audit.csv"
    with audit_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    make_review_sheets(out_rows, args.pic_dir, render_dir, args.out_dir / "review_sheets", args.review_top)
    print(
        {
            "event": "complete",
            "audit": str(audit_path),
            "rows": len(out_rows),
            "rendered": len(list(render_dir.iterdir())) if render_dir.exists() else 0,
            "review_pages": len(list((args.out_dir / "review_sheets").glob("review_page_*.png"))),
        },
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
