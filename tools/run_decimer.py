#!/usr/bin/env python3
import argparse
import csv
import json
import time
from pathlib import Path


def image_files(path: Path):
    return sorted(
        [p for p in path.iterdir() if p.suffix.lower() in {".png", ".jpg", ".jpeg"}],
        key=lambda p: p.name,
    )


def load_existing(path: Path):
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as f:
        return {row["file_name"]: row for row in csv.DictReader(f)}


def wanted_files(files, filter_csv: Path | None, mode: str):
    if filter_csv is None:
        return files
    with filter_csv.open("r", encoding="utf-8", newline="") as f:
        rows = {row["file_name"]: row for row in csv.DictReader(f)}
    selected = []
    for p in files:
        row = rows.get(p.name, {})
        pred = row.get("e_smiles", "")
        try:
            conf = float(row.get("confidence", "0") or 0)
        except ValueError:
            conf = 0.0
        if mode == "all":
            selected.append(p)
        elif mode == "invalid" and (not pred or pred == "<invalid>"):
            selected.append(p)
        elif mode == "weak" and (not pred or pred == "<invalid>" or conf < 0.65):
            selected.append(p)
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--filter-csv", default="")
    parser.add_argument("--mode", choices=["all", "invalid", "weak"], default="all")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    files = image_files(input_dir)
    selected = wanted_files(files, Path(args.filter_csv) if args.filter_csv else None, args.mode)
    if args.limit:
        selected = selected[: args.limit]

    existing = load_existing(out_path) if args.resume else {}
    done = {name for name, row in existing.items() if row.get("smiles") or row.get("error")}
    pending = [p for p in selected if p.name not in done]

    print(json.dumps({"event": "import_decimer", "selected": len(selected), "pending": len(pending)}), flush=True)
    from DECIMER import predict_SMILES

    rows = dict(existing)
    t0 = time.time()
    for idx, p in enumerate(pending, 1):
        tic = time.time()
        pred = ""
        err = ""
        try:
            pred = predict_SMILES(str(p)) or ""
        except Exception as exc:
            err = repr(exc)
        rows[p.name] = {
            "file_name": p.name,
            "smiles": pred,
            "seconds": f"{time.time() - tic:.3f}",
            "error": err,
        }
        with out_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["file_name", "smiles", "seconds", "error"])
            writer.writeheader()
            for q in selected:
                writer.writerow(rows.get(q.name, {"file_name": q.name, "smiles": "", "seconds": "", "error": ""}))
        if idx == 1 or idx % 25 == 0:
            print(
                json.dumps(
                    {
                        "event": "item",
                        "done": idx,
                        "pending": len(pending),
                        "last": p.name,
                        "seconds": round(time.time() - tic, 3),
                        "elapsed": round(time.time() - t0, 3),
                    }
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
