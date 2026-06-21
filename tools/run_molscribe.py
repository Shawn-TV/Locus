#!/usr/bin/env python3
import argparse
import csv
import json
import time
from pathlib import Path

import cv2
import torch
from molscribe import MolScribe
import molscribe.chemistry as chemistry
import molscribe.interface as interface


def convert_graph_to_smiles_seq(coords, symbols, edges, images=None, num_workers=1):
    results = []
    if images is None:
        iterator = zip(coords, symbols, edges)
        for c, s, e in iterator:
            results.append(chemistry._convert_graph_to_smiles(c, s, e))
    else:
        iterator = zip(coords, symbols, edges, images)
        for c, s, e, image in iterator:
            results.append(chemistry._convert_graph_to_smiles(c, s, e, image))
    smiles_list, molblock_list, success = zip(*results)
    return list(smiles_list), list(molblock_list), sum(bool(x) for x in success) / len(success)


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


def predict_files_with_raw(model, paths, batch_size):
    outputs = []
    model.decoder.compute_confidence = True
    for start in range(0, len(paths), batch_size):
        batch_paths = paths[start : start + batch_size]
        images_np = []
        tensors = []
        for path in batch_paths:
            image = cv2.imread(str(path))
            image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
            images_np.append(image)
            tensors.append(model.transform(image=image, keypoints=[])["image"])
        images = torch.stack(tensors, dim=0).to(model.device)
        with torch.no_grad():
            features, hiddens = model.encoder(images)
            predictions = model.decoder.decode(features, hiddens)

        raw_smiles = [pred["chartok_coords"]["smiles"] for pred in predictions]
        node_coords = [pred["chartok_coords"]["coords"] for pred in predictions]
        node_symbols = [pred["chartok_coords"]["symbols"] for pred in predictions]
        edges = [pred["edges"] for pred in predictions]
        post_smiles, molblocks, _ = convert_graph_to_smiles_seq(node_coords, node_symbols, edges, images=images_np)

        for raw, post, molblock, pred in zip(raw_smiles, post_smiles, molblocks, predictions):
            outputs.append(
                {
                    "smiles": post,
                    "raw_smiles": raw,
                    "molfile": molblock,
                    "confidence": pred.get("overall_score", ""),
                }
            )
    return outputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ckpt", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--device", default="cpu", choices=["cpu", "mps"])
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    interface.convert_graph_to_smiles = convert_graph_to_smiles_seq

    input_dir = Path(args.input_dir)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    files = image_files(input_dir)
    if args.limit:
        files = files[: args.limit]

    existing = load_existing(out_path) if args.resume else {}
    done = {name for name, row in existing.items() if row.get("e_smiles")}
    pending = [p for p in files if p.name not in done]

    device = torch.device(args.device)
    print(json.dumps({"event": "start", "device": str(device), "total": len(files), "pending": len(pending)}), flush=True)
    model = MolScribe(args.ckpt, device=device)

    rows = dict(existing)
    t0 = time.time()
    for start in range(0, len(pending), args.batch_size):
        batch = pending[start : start + args.batch_size]
        tic = time.time()
        try:
            outputs = predict_files_with_raw(model, batch, args.batch_size)
        except Exception as exc:
            outputs = [{"smiles": "", "raw_smiles": "", "confidence": "", "error": repr(exc)} for _ in batch]
        for p, output in zip(batch, outputs):
            rows[p.name] = {
                "file_name": p.name,
                "e_smiles": output.get("smiles", "") or "",
                "raw_smiles": output.get("raw_smiles", "") or "",
                "confidence": output.get("confidence", ""),
                "error": output.get("error", ""),
            }
        with out_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["file_name", "e_smiles", "raw_smiles", "confidence", "error"])
            writer.writeheader()
            for p in files:
                writer.writerow(
                    rows.get(
                        p.name,
                        {"file_name": p.name, "e_smiles": "", "raw_smiles": "", "confidence": "", "error": ""},
                    )
                )
        print(
            json.dumps(
                {
                    "event": "batch",
                    "done": min(start + len(batch), len(pending)),
                    "pending": len(pending),
                    "seconds": round(time.time() - tic, 3),
                    "elapsed": round(time.time() - t0, 3),
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
