#!/usr/bin/env python3
"""Batch MolParser web OCSR calls for the Bohrium track.

The public MolParser frontend posts a base64 image to /mol/img2mol and returns
ESMILES captions.  This runner keeps raw response fields so the final combiner
can prefer high-confidence results while preserving a resumable audit trail.
"""

from __future__ import annotations

import argparse
import base64
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable


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


def iter_images(input_dir: Path) -> list[Path]:
    exts = {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff"}
    return sorted(p for p in input_dir.iterdir() if p.suffix.lower() in exts)


def load_done(output: Path) -> set[str]:
    if not output.exists() or output.stat().st_size == 0:
        return set()
    with output.open("r", encoding="utf-8", newline="") as f:
        return {
            row["file_name"]
            for row in csv.DictReader(f)
            if row.get("file_name") and not row.get("error") and (row.get("caption") or row.get("smi"))
        }


def append_rows(output: Path, rows: Iterable[dict[str, object]]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    exists = output.exists() and output.stat().st_size > 0
    with output.open("a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if not exists:
            writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def call_api(path: Path, timeout: float, retries: int, pause: float) -> dict[str, object]:
    started = time.time()
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    payload = json.dumps({"base64_img": encoded}).encode("utf-8")
    headers = {
        "Content-Type": "application/json",
        "User-Agent": "Mozilla/5.0 Bohrium-V1-OCSR",
    }
    last_error = ""

    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(API_URL, data=payload, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                body = response.read().decode("utf-8")
            parsed = json.loads(body)
            if parsed.get("code") != 0:
                raise RuntimeError(parsed.get("msg") or f"nonzero code {parsed.get('code')}")
            results = parsed.get("data") or []
            best = results[0] if results else {}
            return {
                "file_name": path.name,
                "caption": best.get("caption", ""),
                "smi": best.get("smi", ""),
                "score": best.get("score", ""),
                "markush": best.get("markush", ""),
                "n_results": len(results),
                "bbox": json.dumps(best.get("bbox", ""), ensure_ascii=False),
                "trace_id": parsed.get("trace_id", ""),
                "error": "",
                "seconds": round(time.time() - started, 3),
            }
        except (urllib.error.URLError, TimeoutError, RuntimeError, json.JSONDecodeError) as exc:
            last_error = repr(exc)
            if attempt < retries:
                time.sleep(pause * (attempt + 1))

    return {
        "file_name": path.name,
        "caption": "",
        "smi": "",
        "score": "",
        "markush": "",
        "n_results": 0,
        "bbox": "",
        "trace_id": "",
        "error": last_error,
        "seconds": round(time.time() - started, 3),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-pause", type=float, default=2.0)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    images = iter_images(args.input_dir)
    done = load_done(args.output) if args.resume else set()
    pending = [path for path in images if path.name not in done]

    print(
        json.dumps(
            {
                "event": "start",
                "total": len(images),
                "done": len(done),
                "pending": len(pending),
                "workers": args.workers,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )

    completed = 0
    failed = 0
    started = time.time()
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {
            pool.submit(call_api, path, args.timeout, args.retries, args.retry_pause): path
            for path in pending
        }
        buffer: list[dict[str, object]] = []
        for future in as_completed(futures):
            row = future.result()
            buffer.append(row)
            completed += 1
            failed += 1 if row.get("error") else 0
            if len(buffer) >= 20:
                append_rows(args.output, buffer)
                buffer.clear()
            if completed % 20 == 0 or completed == len(pending):
                print(
                    json.dumps(
                        {
                            "event": "progress",
                            "completed": completed,
                            "pending": len(pending),
                            "failed": failed,
                            "elapsed": round(time.time() - started, 3),
                        },
                        ensure_ascii=False,
                    ),
                    flush=True,
                )
        if buffer:
            append_rows(args.output, buffer)

    print(
        json.dumps(
            {
                "event": "complete",
                "completed": completed,
                "failed": failed,
                "elapsed": round(time.time() - started, 3),
            },
            ensure_ascii=False,
        ),
        flush=True,
    )
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
