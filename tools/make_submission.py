#!/usr/bin/env python3
"""Create a competition-ready V1 zip from OCSR model outputs."""

from __future__ import annotations

import argparse
import csv
import re
import zipfile
from collections import Counter
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except Exception:  # pragma: no cover - generation still works without RDKit.
    Chem = None


EXPAND_LABELS = {
    "Me": "C",
    "Et": "CC",
    "Pr": "CCC",
    "nPr": "CCC",
    "iPr": "C(C)C",
    "Bu": "CCCC",
    "nBu": "CCCC",
    "tBu": "C(C)(C)C",
    "Ph": "c1ccccc1",
    "Bn": "Cc1ccccc1",
    "CF3": "C(F)(F)F",
    "OMe": "OC",
    "MeO": "OC",
    "OEt": "OCC",
    "EtO": "OCC",
    "NMe2": "N(C)C",
    "NEt2": "N(CC)CC",
    "NO2": "[N+](=O)[O-]",
    "CN": "C#N",
    "CO2Me": "C(=O)OC",
    "COOMe": "C(=O)OC",
    "CO2Et": "C(=O)OCC",
    "COOEt": "C(=O)OCC",
}

ANNOTATE_LABELS = {
    "R",
    "R1",
    "R2",
    "R3",
    "R4",
    "R5",
    "Ar",
    "Tf",
    "TMS",
    "TBS",
    "TBDMS",
    "TBDPS",
    "Boc",
    "Cbz",
    "Alloc",
    "Fmoc",
    "Ts",
    "Ms",
    "Ns",
    "TIPS",
    "SEM",
    "MEM",
}

ELEMENTS = {
    "B",
    "C",
    "N",
    "O",
    "P",
    "S",
    "F",
    "I",
    "H",
    "Cl",
    "Br",
    "Si",
    "Na",
    "Li",
    "K",
    "Mg",
    "Ca",
    "Zn",
    "Fe",
    "Cu",
    "Al",
    "Se",
    "Sn",
    "Hg",
    "Ag",
    "Au",
    "Pt",
    "Pd",
    "Co",
    "Ni",
}


def read_csv(path: Path, encoding: str = "utf-8") -> list[dict[str, str]]:
    if not path or not path.exists():
        return []
    with path.open("r", encoding=encoding, newline="") as f:
        return list(csv.DictReader(f))


def ensure_esmi(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return ""
    value = re.sub(r"<r>(.*?)</r>", r"<a>\1</a>", value)
    return value if "<sep>" in value else f"{value}<sep>"


def valid_base_score(esmi: str) -> int:
    if not esmi:
        return 0
    if Chem is None:
        return 1
    base = esmi.split("<sep>")[0]
    try:
        return 1 if Chem.MolFromSmiles(base) is not None else 0
    except Exception:
        return 0


def candidate_rank(esmi: str, score: float) -> tuple:
    valid = valid_base_score(esmi)
    if score < 0.5:
        return (valid, 0, len(esmi), score)
    if score < 0.95:
        return (valid, 1, score, len(esmi))
    return (valid, 2, score, len(esmi))


def looks_like_atom_bracket(token: str) -> bool:
    token = token.strip()
    if not token:
        return False
    if token == "*":
        return True
    atom = re.match(r"^\d*([A-Z][a-z]?|[cnopsb]|\*)", token)
    if not atom or (atom.group(1) not in ELEMENTS and atom.group(1) not in {"c", "n", "o", "p", "s", "b", "*"}):
        return False
    rest = token[atom.end() :]
    return bool(re.fullmatch(r"(@@?|H\d?|\d+|[+-]\d*|:\d+|%\d+)*", rest))


def count_atoms(smiles: str) -> int:
    count = 0
    i = 0
    while i < len(smiles):
        ch = smiles[i]
        if ch == "[":
            end = smiles.find("]", i + 1)
            if end == -1:
                i += 1
                continue
            if looks_like_atom_bracket(smiles[i + 1 : end]):
                count += 1
            i = end + 1
            continue
        if ch == "*":
            count += 1
            i += 1
            continue
        if smiles.startswith(("Cl", "Br", "Si", "Na", "Li"), i):
            count += 1
            i += 2
            continue
        if ch in "BCNOPSFIK" or ch in "cnopsb":
            count += 1
        i += 1
    return count


def raw_to_esmi(raw: str) -> str:
    raw = (raw or "").strip().replace(" ", "")
    if not raw:
        return ""

    output: list[str] = []
    annotations: list[tuple[int, str]] = []
    atom_index = 0
    i = 0
    while i < len(raw):
        if raw[i] == "[":
            end = raw.find("]", i + 1)
            if end != -1:
                token = raw[i + 1 : end]
                if token in EXPAND_LABELS:
                    replacement = EXPAND_LABELS[token]
                    output.append(replacement)
                    atom_index += count_atoms(replacement)
                    i = end + 1
                    continue
                if token in ANNOTATE_LABELS or re.match(r"^(R|R\d+|[A-Z][A-Za-z0-9]*\[[0-9]+\])$", token):
                    output.append("*")
                    annotations.append((atom_index, token))
                    atom_index += 1
                    i = end + 1
                    continue
                if looks_like_atom_bracket(token):
                    output.append(raw[i : end + 1])
                    atom_index += 1
                    i = end + 1
                    continue
                output.append("*")
                annotations.append((atom_index, token))
                atom_index += 1
                i = end + 1
                continue
        if raw.startswith(("Cl", "Br", "Si", "Na", "Li"), i):
            output.append(raw[i : i + 2])
            atom_index += 1
            i += 2
            continue
        output.append(raw[i])
        if raw[i] == "*" or raw[i] in "BCNOPSFIKcnopsb":
            atom_index += 1
        i += 1

    base = "".join(output)
    suffix = "".join(f"<a>{idx}:{label}</a>" for idx, label in annotations)
    return f"{base}<sep>{suffix}"


def load_molparser(paths: Path | list[Path]) -> dict[str, tuple[str, float]]:
    if isinstance(paths, Path):
        paths = [paths]
    best: dict[str, tuple[str, float, tuple[int, float, int]]] = {}
    for path in paths:
        for row in read_csv(path):
            if row.get("error"):
                continue
            file_name = row.get("file_name", "")
            candidate = ensure_esmi(row.get("caption", "")) or ensure_esmi(row.get("smi", ""))
            if not file_name or not candidate:
                continue
            try:
                score = float(row.get("score") or 0.0)
            except ValueError:
                score = 0.0
            rank = candidate_rank(candidate, score)
            previous = best.get(file_name)
            if previous is None or rank > previous[2]:
                best[file_name] = (candidate, score, rank)
    return {file_name: (candidate, score) for file_name, (candidate, score, _rank) in best.items()}


def load_molscribe(path: Path) -> dict[str, tuple[str, str]]:
    outputs: dict[str, tuple[str, str]] = {}
    for row in read_csv(path):
        if row.get("error"):
            continue
        file_name = row.get("file_name", "")
        if not file_name:
            continue
        post = ensure_esmi(row.get("e_smiles", ""))
        raw = raw_to_esmi(row.get("raw_smiles", ""))
        if post and "<invalid>" not in post:
            outputs[file_name] = (post, "molscribe")
        elif raw:
            outputs[file_name] = (raw, "molscribe_raw")
    return outputs


def load_overrides(path: Path | None) -> dict[str, str]:
    if not path:
        return {}
    rows = read_csv(path, encoding="utf-8-sig")
    return {
        row["file_name"]: ensure_esmi(row["e_smiles"])
        for row in rows
        if row.get("file_name") and row.get("e_smiles")
    }


def template_names(template_csv: Path, pic_dir: Path) -> list[str]:
    rows = read_csv(template_csv, encoding="utf-8-sig")
    if rows and "file_name" in rows[0]:
        return [row["file_name"] for row in rows if row.get("file_name")]
    exts = {".png", ".jpg", ".jpeg"}
    return sorted(p.name for p in pic_dir.iterdir() if p.suffix.lower() in exts)


def write_meta(
    path: Path,
    stats: Counter[str],
    molparser_count: int,
    molscribe_count: int,
    version: str,
    code_repo: str,
) -> None:
    total = sum(stats.values())
    text = f"""# Track 1 submission notes (meta.md)

## 1. Model and compute usage

| Model / tool | Source | Version | Invocation | Token usage | Notes |
|--------------|--------|---------|------------|-------------|-------|
| MolParser | DP Technology | Public OCSR service | API / Web service | N/A | Primary OCSR candidate source for original and preprocessed images |
| MolScribe | Open-source checkpoint | User-provided | Local inference | N/A | Optional fallback for missing or difficult rows |
| RDKit | Open-source | Local Python package | Local validation/rendering | N/A | SMILES validation and rendered review sheets |

**Total model token usage**: N/A

## 2. External data and code declaration

- **Code repository**: {code_repo}

External models / tools used by the workflow:

- MolParser public OCSR service by DP Technology.
- Optional MolScribe released checkpoint, subject to the upstream project license.
- RDKit, subject to the upstream RDKit license.
- DECIMER runner is included as an optional helper; use only if installed separately under its upstream license.

This repository does not include the competition test images, private submission CSV files, model weights, or proprietary third-party service code.

## 3. Annotation / prediction method

Version: {version}

The submission is generated for Bohrium competition 53859761357, Track 1.

Pipeline summary:

- Primary recognition: MolParser public OCSR service, preserving returned E-SMILES captions.
- Image preprocessing: original image plus selected crop/pad/upscale, binary, and color-foreground cleanup variants for weak recognitions.
- Candidate selection: prefer valid underlying SMILES, high-confidence MolParser results, and fuller low-confidence structures.
- Optional fallback: local MolScribe inference with light E-SMILES normalization for missing online results.
- Audit: RDKit renders predictions and produces review sheets for human visual checking.
- Formatting: UTF-8 CSV with required columns `file_name,e_smiles`; all populated rows end with `<sep>` or include E-SMILES extension tags.

Coverage summary:

- Total rows: {total}
- Primary MolParser rows used: {stats.get('molparser', 0)}
- Manual visual correction rows used: {stats.get('manual', 0)}
- MolScribe fallback rows used: {stats.get('molscribe', 0) + stats.get('molscribe_raw', 0)}
- Empty rows: {stats.get('empty', 0)}
- MolParser unique successful recognitions available: {molparser_count}
- MolScribe unique recognitions available: {molscribe_count}
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pic-dir", required=True, type=Path)
    parser.add_argument("--template-csv", required=True, type=Path)
    parser.add_argument("--molparser-csv", required=True, nargs="+", type=Path)
    parser.add_argument("--molscribe-csv", required=True, type=Path)
    parser.add_argument("--overrides-csv", type=Path)
    parser.add_argument("--version", default="V1")
    parser.add_argument("--code-repo", default="N/A")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--zip-path", required=True, type=Path)
    args = parser.parse_args()

    names = template_names(args.template_csv, args.pic_dir)
    molparser = load_molparser(args.molparser_csv)
    molscribe = load_molscribe(args.molscribe_csv)
    overrides = load_overrides(args.overrides_csv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    submission_path = args.out_dir / "submission.csv"
    meta_path = args.out_dir / "meta.md"

    stats: Counter[str] = Counter()
    with submission_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["file_name", "e_smiles"])
        writer.writeheader()
        for name in names:
            if name in overrides:
                value = overrides[name]
                source = "manual"
            elif name in molparser:
                value = molparser[name][0]
                source = "molparser"
            elif name in molscribe:
                value, source = molscribe[name]
            else:
                value = ""
                source = "empty"
            stats[source] += 1
            writer.writerow({"file_name": name, "e_smiles": value})

    write_meta(meta_path, stats, len(molparser), len(molscribe), args.version, args.code_repo)

    args.zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(submission_path, "submission.csv")
        zf.write(meta_path, "meta.md")

    print(
        {
            "zip": str(args.zip_path),
            "rows": len(names),
            "stats": dict(stats),
            "molparser_success": len(molparser),
            "molscribe_success": len(molscribe),
        }
    )
    return 0 if stats.get("empty", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
