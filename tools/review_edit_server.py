#!/usr/bin/env python3
"""Local review/edit UI for E-SMILES predictions.

The server is intentionally local-only by default. It loads a base submission,
lets a reviewer edit one molecule at a time, renders the edited E-SMILES with
RDKit, and saves changes to both a full working CSV and a compact overrides CSV.
"""

from __future__ import annotations

import argparse
import csv
import html
import io
import json
import os
import re
import sys
import threading
import urllib.parse
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

if os.name == "nt" and hasattr(os, "add_dll_directory"):
    site_packages = Path(sys.executable).parent / "Lib" / "site-packages"
    for dll_dir in site_packages.glob("*.libs"):
        if dll_dir.is_dir():
            os.add_dll_directory(str(dll_dir))

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


RDLogger.DisableLog("rdApp.*")
ANNOT_RE = re.compile(r"<([ar])>(\d+):(.*?)</\1>")


def translate_candidate_source(source: str) -> str:
    labels = {
        "molparser_molparser_variants_safe": "MolParser 安全预处理",
        "molparser_molparser_variants_aggressive": "MolParser 强预处理",
        "molparser_molparser_variants_color_low": "MolParser 降低彩色干扰",
        "molparser_molparser_variants_v5_nonmp": "非 MolParser 候选",
    }
    return labels.get(source, source)


def translate_reason_piece(reason: str) -> str:
    reason = reason.strip()
    if not reason:
        return ""
    match = re.fullmatch(r"low original score ([0-9.]+)", reason)
    if match:
        return f"原始 MolParser 置信度低（{match.group(1)}）"
    match = re.fullmatch(r"medium-low original score ([0-9.]+)", reason)
    if match:
        return f"原始 MolParser 置信度中低（{match.group(1)}）"
    match = re.fullmatch(r"weak original score ([0-9.]+)", reason)
    if match:
        return f"原始 MolParser 置信度偏弱（{match.group(1)}）"
    match = re.fullmatch(r"sub-0.95 original score ([0-9.]+)", reason)
    if match:
        return f"原始 MolParser 分数低于 0.95（{match.group(1)}）"
    if reason == "no original MolParser success":
        return "原始 MolParser 没有成功返回结果"
    match = re.fullmatch(r"high-score different candidate ([^:]+):([0-9.]+)", reason)
    if match:
        return f"存在高分但不同的候选结果（{translate_candidate_source(match.group(1))}，分数 {match.group(2)}）"
    match = re.fullmatch(r"different candidate ([^:]+):([0-9.]+)", reason)
    if match:
        return f"存在另一个不同候选结果（{translate_candidate_source(match.group(1))}，分数 {match.group(2)}）"
    if reason == "has E-SMILES labels":
        return "原图包含 R、Ph、Me 等 E-SMILES 标注，容易出现基团或连接位识别错误"
    match = re.fullmatch(r"colored/label-heavy image ([0-9.]+)", reason)
    if match:
        return f"原图彩色信息或文字标签较多，可能干扰识别（指标 {match.group(1)}）"
    match = re.fullmatch(r"render visual mismatch phash ([0-9.]+)", reason)
    if match:
        return f"当前结果重画后与原图视觉差异较大（pHash 差异 {match.group(1)}）"
    match = re.fullmatch(r"ink density delta ([0-9.]+)", reason)
    if match:
        return f"原图与重画图线条密度差异较大（差值 {match.group(1)}）"
    match = re.fullmatch(r"low render SSIM ([0-9.]+)", reason)
    if match:
        return f"重画图与原图结构相似度较低（SSIM {match.group(1)}）"
    if reason == "very small/label-like structure":
        return "结构很小或接近标签图，容易被识别成文字/基团而不是完整结构"
    return reason


def translate_reasons(reasons: str, is_priority: bool, priority_count: int = 300) -> str:
    if not reasons.strip():
        if is_priority:
            return "没有记录到具体风险原因。"
        if priority_count <= 0:
            return "没有加载风险队列；这张仍在全量工作稿里，可以搜索、修改，并会随工作 zip 一起导出。"
        return f"未进入前 {priority_count} 重点风险队列；这张仍在全量工作稿里，可以搜索、修改，并会随工作 zip 一起导出。"
    translated = [translate_reason_piece(part) for part in reasons.split(";")]
    translated = [part for part in translated if part]
    return "\n".join(f"{idx}. {text}" for idx, text in enumerate(translated, 1))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def ensure_esmi(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return ""
    return value if "<sep>" in value else f"{value}<sep>"


def base_part(esmi: str) -> str:
    return ensure_esmi(esmi).split("<sep>")[0]


def annotation_items(esmi: str, kind: str) -> list[tuple[int, str]]:
    out: list[tuple[int, str]] = []
    for annot_kind, idx, label in ANNOT_RE.findall(ensure_esmi(esmi)):
        if annot_kind == kind:
            out.append((int(idx), label))
    return out


def atom_label_for_draw(label: str) -> str:
    return "*" if label == "<dum>" else label


def render_error(message: str, width: int = 780, height: int = 560) -> bytes:
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.rectangle((12, 12, width - 13, height - 13), outline=(220, 70, 70), width=3)
    draw.text((32, 32), "Render failed", fill=(160, 0, 0))
    y = 70
    for line in message[:900].splitlines() or [message[:900]]:
        while len(line) > 82:
            draw.text((32, y), line[:82], fill=(50, 50, 50))
            line = line[82:]
            y += 22
        draw.text((32, y), line, fill=(50, 50, 50))
        y += 24
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def render_esmi_png(esmi: str, size: tuple[int, int] = (900, 640)) -> bytes:
    esmi = ensure_esmi(esmi)
    base = base_part(esmi)
    if not base:
        return render_error("Empty E-SMILES")
    mol = Chem.MolFromSmiles(base)
    if mol is None:
        return render_error(f"RDKit cannot parse base SMILES:\n{base}")
    mol = Chem.Mol(mol)
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    options = drawer.drawOptions()
    options.clearBackground = True
    options.bondLineWidth = 2.2
    atom_labels: dict[int, list[str]] = {}
    for idx, label in annotation_items(esmi, "a"):
        if 0 <= idx < mol.GetNumAtoms():
            atom_labels.setdefault(idx, []).append(atom_label_for_draw(label))
    for idx, labels in atom_labels.items():
        options.atomLabels[idx] = "/".join(labels)

    try:
        drawer.DrawMolecule(mol)
    except Exception as exc:
        return render_error(f"RDKit drawing error:\n{exc}\n\nBase SMILES:\n{base}")

    ring_labels: dict[int, list[str]] = {}
    for ring_idx, label in annotation_items(esmi, "r"):
        ring_labels.setdefault(ring_idx, []).append(label)

    ring_label_positions: list[tuple[float, float, list[str]]] = []
    rings = mol.GetRingInfo().AtomRings()
    for ring_idx, labels in ring_labels.items():
        if 0 <= ring_idx < len(rings) and rings[ring_idx]:
            points = [drawer.GetDrawCoords(atom_idx) for atom_idx in rings[ring_idx]]
            x = sum(point.x for point in points) / len(points)
            y = sum(point.y for point in points) / len(points)
            ring_label_positions.append((x, y, labels))

    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    if not ring_label_positions:
        return png

    img = Image.open(io.BytesIO(png)).convert("RGB")
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 24)
    except Exception:
        font = None
    for x, y, labels in ring_label_positions:
        line_metrics: list[tuple[str, int, int]] = []
        for label in labels:
            bbox = draw.textbbox((0, 0), label, font=font)
            line_metrics.append((label, bbox[2] - bbox[0], bbox[3] - bbox[1]))
        line_height = max([height for _, _, height in line_metrics] + [18]) + 4
        start_y = y - (line_height * len(line_metrics)) / 2
        for offset, (label, width, _height) in enumerate(line_metrics):
            draw.text((x - width / 2, start_y + offset * line_height), label, fill=(185, 0, 0), font=font)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def image_preview_png(src: Path, cache: Path, max_side: int = 1600) -> bytes:
    if cache.exists() and cache.stat().st_mtime >= src.stat().st_mtime:
        return cache.read_bytes()
    cache.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(src) as img:
        img = img.copy()
        if max(img.size) > max_side:
            try:
                resampling = Image.Resampling.LANCZOS
            except AttributeError:
                resampling = Image.LANCZOS
            img.thumbnail((max_side, max_side), resampling)
        if img.mode not in ("RGB", "RGBA"):
            img = img.convert("RGB")
        buf = io.BytesIO()
        img.save(buf, format="PNG", optimize=True)
    data = buf.getvalue()
    cache.write_bytes(data)
    return data


class ReviewState:
    def __init__(self, args: argparse.Namespace) -> None:
        self.pic_dir: Path = args.pic_dir
        self.output_dir: Path = args.output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.working_csv = self.output_dir / "working_submission.csv"
        self.overrides_csv = self.output_dir / "review_overrides.csv"
        self.review_state_json = self.output_dir / "review_state.json"
        self.preview_dir = self.output_dir / "image_preview_cache"
        self.package_zip = self.output_dir / "working_submission_package.zip"
        self.lock = threading.Lock()

        self.rows = read_csv(args.submission_csv)
        self.order = [row["file_name"] for row in self.rows]
        self.values = {row["file_name"]: ensure_esmi(row["e_smiles"]) for row in self.rows}
        self.original_values = dict(self.values)

        if self.working_csv.exists():
            for row in read_csv(self.working_csv):
                if row.get("file_name") in self.values:
                    self.values[row["file_name"]] = ensure_esmi(row.get("e_smiles", ""))

        self.risk_rows = read_csv(args.risk_csv) if args.risk_csv and args.risk_csv.exists() else []
        self.priority_count = len(self.risk_rows)
        known = {row.get("file_name", "") for row in self.risk_rows}
        for name in self.order:
            if name not in known:
                self.risk_rows.append(
                    {
                        "rank": "",
                        "page": "",
                        "slot": "",
                        "file_name": name,
                        "risk": "",
                        "reasons": "",
                        "orig_score": "",
                        "best_source": "",
                        "best_score": "",
                        "diff_source": "",
                        "diff_score": "",
                        "v1_atoms": "",
                        "v1_esmi": self.original_values.get(name, ""),
                        "diff_candidate": "",
                    }
                )

        self.meta_path: Path | None = args.meta_path if args.meta_path and args.meta_path.exists() else None
        self.save_files()

    def read_review_state(self) -> dict[str, Any]:
        if not self.review_state_json.exists():
            return {}
        try:
            data = json.loads(self.review_state_json.read_text(encoding="utf-8"))
        except Exception:
            return {}
        return data if isinstance(data, dict) else {}

    def save_review_state(self, data: dict[str, Any]) -> None:
        clean: dict[str, dict[str, Any]] = {}
        for name, value in data.items():
            if not isinstance(name, str) or name not in self.values or not isinstance(value, dict):
                continue
            clean[name] = {
                "reviewed": bool(value.get("reviewed")),
                "flagged": bool(value.get("flagged")),
                "note": str(value.get("note", "")),
            }
        tmp = self.review_state_json.with_suffix(".tmp")
        tmp.write_text(json.dumps(clean, ensure_ascii=False, indent=2), encoding="utf-8")
        tmp.replace(self.review_state_json)

    def item_payload(self, row: dict[str, str]) -> dict[str, str]:
        name = row["file_name"]
        payload = dict(row)
        rank_text = (row.get("rank") or "").strip()
        is_priority = rank_text.isdigit() and int(rank_text) <= self.priority_count
        payload["is_priority_review"] = str(is_priority).lower()
        payload["priority_label"] = "重点复查" if is_priority else "全量样本"
        if is_priority:
            payload["priority_note"] = f"前 {self.priority_count} 张重点风险队列，建议优先人工核查。"
        elif self.priority_count:
            payload["priority_note"] = f"未进入前 {self.priority_count} 张重点风险队列；仍包含在全量导出中。"
        else:
            payload["priority_note"] = "没有加载风险队列；仍包含在全量导出中。"
        payload["risk_total"] = str(self.priority_count)
        payload["reasons_zh"] = translate_reasons(row.get("reasons", ""), is_priority, self.priority_count)
        payload["best_source_zh"] = translate_candidate_source(row.get("best_source", ""))
        payload["diff_source_zh"] = translate_candidate_source(row.get("diff_source", ""))
        payload["e_smiles"] = self.values.get(name, "")
        payload["original_e_smiles"] = self.original_values.get(name, "")
        payload["changed"] = str(payload["e_smiles"] != payload["original_e_smiles"]).lower()
        return payload

    def save_files(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        with self.working_csv.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["file_name", "e_smiles"])
            writer.writeheader()
            for name in self.order:
                writer.writerow({"file_name": name, "e_smiles": self.values.get(name, "")})

        with self.overrides_csv.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["file_name", "e_smiles", "original_e_smiles"])
            writer.writeheader()
            for name in self.order:
                value = self.values.get(name, "")
                original = self.original_values.get(name, "")
                if value != original:
                    writer.writerow({"file_name": name, "e_smiles": value, "original_e_smiles": original})

    def build_zip(self) -> None:
        meta_text = self.meta_path.read_text(encoding="utf-8") if self.meta_path else "# meta.md\n"
        with zipfile.ZipFile(self.package_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.write(self.working_csv, "submission.csv")
            zf.writestr("meta.md", meta_text)


INDEX_HTML = r"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Locus Review Editor</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f4f6f8;
      --line: #d8dee8;
      --text: #131821;
      --muted: #5d6676;
      --blue: #1456c8;
      --green: #137333;
      --red: #b3261e;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 14px;
    }
    header {
      height: 54px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      padding: 0 16px;
      background: white;
      border-bottom: 1px solid var(--line);
    }
    h1 { margin: 0; font-size: 17px; }
    .status { color: var(--muted); font-size: 13px; }
    .app {
      height: calc(100vh - 54px);
      display: grid;
      grid-template-columns: 330px minmax(0, 1fr);
      min-height: 680px;
    }
    aside {
      border-right: 1px solid var(--line);
      background: #fbfcfe;
      overflow: hidden;
      display: flex;
      flex-direction: column;
    }
    .filters {
      padding: 10px;
      display: grid;
      grid-template-columns: 1fr auto auto;
      gap: 8px;
      border-bottom: 1px solid var(--line);
    }
    input, textarea, button {
      font: inherit;
    }
    input, textarea {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: white;
      color: var(--text);
    }
    input {
      height: 34px;
      padding: 0 9px;
      min-width: 0;
    }
    button {
      height: 34px;
      padding: 0 11px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: white;
      cursor: pointer;
    }
    button.primary {
      background: var(--blue);
      border-color: var(--blue);
      color: white;
      font-weight: 650;
    }
    button:disabled {
      cursor: default;
      opacity: 0.55;
    }
    .list {
      overflow: auto;
      padding: 6px;
    }
    .item {
      width: 100%;
      height: auto;
      display: block;
      text-align: left;
      padding: 8px;
      margin: 0 0 5px;
      border-radius: 6px;
      background: white;
    }
    .item.active {
      border-color: var(--blue);
      outline: 2px solid rgba(20, 86, 200, 0.18);
    }
    .item.changed::after {
      content: "saved";
      float: right;
      color: var(--green);
      font-size: 11px;
      font-weight: 700;
    }
    .item-title {
      display: flex;
      justify-content: space-between;
      gap: 8px;
      font-weight: 700;
    }
    .item-meta {
      margin-top: 3px;
      color: var(--muted);
      font-size: 12px;
      line-height: 1.35;
    }
    main {
      min-width: 0;
      overflow: auto;
      padding: 12px;
    }
    .topbar {
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
      margin-bottom: 10px;
    }
    .title {
      font-size: 16px;
      font-weight: 750;
      margin-right: auto;
    }
    .compare {
      display: grid;
      grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
      gap: 12px;
      min-height: 520px;
    }
    .pane {
      background: white;
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      display: flex;
      flex-direction: column;
      min-width: 0;
    }
    .pane h2 {
      margin: 0;
      padding: 9px 11px;
      font-size: 13px;
      border-bottom: 1px solid var(--line);
      background: #fbfcfe;
    }
    .image-wrap {
      flex: 1;
      min-height: 430px;
      display: flex;
      align-items: center;
      justify-content: center;
      background: white;
      overflow: auto;
      padding: 12px;
    }
    .image-wrap img {
      max-width: 100%;
      max-height: 76vh;
      object-fit: contain;
    }
    .editor {
      margin-top: 12px;
      display: grid;
      grid-template-columns: minmax(0, 1fr) 320px;
      gap: 12px;
    }
    textarea {
      width: 100%;
      min-height: 112px;
      resize: vertical;
      padding: 10px;
      line-height: 1.45;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 13px;
    }
    .info {
      background: white;
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 10px;
      color: var(--muted);
      overflow-wrap: anywhere;
      line-height: 1.45;
    }
    .info strong { color: var(--text); }
    .save-state.ok { color: var(--green); }
    .save-state.dirty { color: var(--red); }
    .candidate {
      margin-top: 8px;
      padding: 8px;
      background: #f8fafd;
      border: 1px solid #e4e9f2;
      border-radius: 6px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 12px;
      color: #334155;
    }
    @media (max-width: 1100px) {
      .app { grid-template-columns: 260px minmax(0, 1fr); }
      .compare, .editor { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
  <header>
    <h1>Locus Review Editor</h1>
    <div class="status" id="globalStatus">loading</div>
  </header>
  <div class="app">
    <aside>
      <div class="filters">
        <input id="search" placeholder="file / rank / reason">
        <button id="changedOnly">Changed</button>
        <button id="allItems">All</button>
      </div>
      <div class="list" id="list"></div>
    </aside>
    <main>
      <div class="topbar">
        <div class="title" id="title">No item selected</div>
        <button id="prevBtn">Prev</button>
        <button id="nextBtn">Next</button>
        <button id="resetBtn">Reset</button>
        <button class="primary" id="saveBtn">Save</button>
        <button id="zipBtn">Build zip</button>
        <span class="save-state" id="saveState"></span>
      </div>
      <div class="compare">
        <section class="pane">
          <h2>Original image</h2>
          <div class="image-wrap"><img id="origImg" alt="original"></div>
        </section>
        <section class="pane">
          <h2>Rendered from edited E-SMILES</h2>
          <div class="image-wrap"><img id="renderImg" alt="rendered"></div>
        </section>
      </div>
      <div class="editor">
        <textarea id="editor" spellcheck="false"></textarea>
        <div class="info">
          <div><strong id="fileName">-</strong></div>
          <div id="riskInfo"></div>
          <div id="reasonInfo"></div>
          <div id="scoreInfo"></div>
          <div class="candidate" id="candidateInfo"></div>
        </div>
      </div>
    </main>
  </div>
  <script>
    let items = [];
    let filtered = [];
    let currentIndex = 0;
    let current = null;
    let dirty = false;
    let changedOnly = false;
    let renderTimer = null;

    const $ = (id) => document.getElementById(id);

    function setStatus(text) {
      $('globalStatus').textContent = text;
    }

    function pageSlot(item) {
      const page = item.page || '-';
      const slot = item.slot || '-';
      return `page ${page}-${slot}`;
    }

    function applyFilter() {
      const q = $('search').value.trim().toLowerCase();
      filtered = items.filter((item) => {
        if (changedOnly && item.changed !== 'true') return false;
        if (!q) return true;
        return [item.file_name, item.rank, item.reasons, item.e_smiles, item.diff_candidate]
          .join(' ').toLowerCase().includes(q);
      });
      renderList();
      if (filtered.length && (!current || !filtered.find((x) => x.file_name === current.file_name))) {
        selectByIndex(0);
      }
    }

    function renderList() {
      const list = $('list');
      list.innerHTML = '';
      filtered.forEach((item, index) => {
        const btn = document.createElement('button');
        btn.className = 'item';
        if (current && item.file_name === current.file_name) btn.classList.add('active');
        if (item.changed === 'true') btn.classList.add('changed');
        btn.innerHTML = `
          <div class="item-title"><span>#${item.rank || '?' } ${item.file_name}</span><span>${Number(item.risk || 0).toFixed(1)}</span></div>
          <div class="item-meta">${pageSlot(item)} · score ${item.orig_score || '-'}</div>
          <div class="item-meta">${escapeHtml((item.reasons || '').slice(0, 110))}</div>
        `;
        btn.addEventListener('click', () => selectByIndex(index));
        list.appendChild(btn);
      });
    }

    function escapeHtml(s) {
      return s.replace(/[&<>"']/g, (c) => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
    }

    async function selectByIndex(index) {
      if (index < 0 || index >= filtered.length) return;
      currentIndex = index;
      current = filtered[index];
      dirty = false;
      $('title').textContent = `#${current.rank || '?'} · ${current.file_name} · ${pageSlot(current)}`;
      $('fileName').textContent = current.file_name;
      $('riskInfo').innerHTML = `<strong>Risk:</strong> ${current.risk || '-'} · <strong>Rank:</strong> ${current.rank || '-'}`;
      $('reasonInfo').innerHTML = `<strong>Reasons:</strong> ${escapeHtml(current.reasons || '-')}`;
      $('scoreInfo').innerHTML = `<strong>Original score:</strong> ${current.orig_score || '-'} · <strong>Best:</strong> ${current.best_source || '-'} ${current.best_score || ''}`;
      $('candidateInfo').textContent = current.diff_candidate ? `Different candidate:\n${current.diff_candidate}` : 'No different candidate recorded.';
      $('editor').value = current.e_smiles || '';
      $('origImg').src = `/original/${encodeURIComponent(current.file_name)}`;
      updateSaveState();
      renderCurrent();
      renderList();
    }

    function updateSaveState(text) {
      const state = $('saveState');
      if (text) {
        state.textContent = text;
        state.className = text.includes('Saved') || text.includes('Built') ? 'save-state ok' : 'save-state dirty';
        return;
      }
      state.textContent = dirty ? 'Unsaved changes' : 'Saved';
      state.className = dirty ? 'save-state dirty' : 'save-state ok';
    }

    function scheduleRender() {
      dirty = true;
      updateSaveState();
      clearTimeout(renderTimer);
      renderTimer = setTimeout(renderCurrent, 350);
    }

    async function renderCurrent() {
      if (!current) return;
      const esmi = $('editor').value;
      const res = await fetch('/api/render', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({e_smiles: esmi})
      });
      const blob = await res.blob();
      const old = $('renderImg').src;
      $('renderImg').src = URL.createObjectURL(blob);
      if (old.startsWith('blob:')) URL.revokeObjectURL(old);
    }

    async function saveCurrent() {
      if (!current) return;
      const esmi = $('editor').value.trim();
      const res = await fetch('/api/save', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({file_name: current.file_name, e_smiles: esmi})
      });
      const data = await res.json();
      if (!data.ok) {
        updateSaveState(`Save failed: ${data.error || 'unknown error'}`);
        return;
      }
      current.e_smiles = data.e_smiles;
      current.changed = String(data.changed);
      const global = items.find((x) => x.file_name === current.file_name);
      if (global) Object.assign(global, current);
      dirty = false;
      updateSaveState(`Saved ${new Date().toLocaleTimeString()}`);
      renderList();
    }

    async function buildZip() {
      const res = await fetch('/api/build_zip', {method: 'POST'});
      const data = await res.json();
      updateSaveState(data.ok ? `Built zip: ${data.path}` : `Zip failed: ${data.error}`);
    }

    async function loadItems() {
      const res = await fetch('/api/items');
      items = await res.json();
      filtered = items.slice();
      setStatus(`${items.length} review items loaded`);
      renderList();
      selectByIndex(0);
    }

    $('editor').addEventListener('input', scheduleRender);
    $('saveBtn').addEventListener('click', saveCurrent);
    $('zipBtn').addEventListener('click', buildZip);
    $('prevBtn').addEventListener('click', () => selectByIndex(currentIndex - 1));
    $('nextBtn').addEventListener('click', () => selectByIndex(currentIndex + 1));
    $('resetBtn').addEventListener('click', () => {
      if (!current) return;
      $('editor').value = current.original_e_smiles || '';
      scheduleRender();
    });
    $('search').addEventListener('input', applyFilter);
    $('changedOnly').addEventListener('click', () => { changedOnly = true; applyFilter(); });
    $('allItems').addEventListener('click', () => { changedOnly = false; $('search').value = ''; applyFilter(); });
    window.addEventListener('keydown', (event) => {
      if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === 's') {
        event.preventDefault();
        saveCurrent();
      }
    });
    loadItems().catch((err) => setStatus(`failed: ${err}`));
  </script>
</body>
</html>
"""

WORKBENCH_HTML = r"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Locus 人工复查工作台</title>
  <style>
    :root {
      --bg: #f5f3ef;
      --rail: #ece8df;
      --panel: #fffdf8;
      --ink: #172026;
      --muted: #66717a;
      --line: #d8d5cc;
      --accent: #176b87;
      --accent-2: #2f8f6b;
      --warm: #d98a2b;
      --danger: #b8403a;
      --shadow: 0 16px 50px rgba(31, 39, 45, 0.12);
      --radius: 8px;
      color-scheme: light;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", "PingFang SC", "Microsoft YaHei", sans-serif;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      min-width: 1180px;
      height: 100vh;
      overflow: hidden;
      background: var(--bg);
      color: var(--ink);
      font-size: 14px;
    }
    button, input, textarea { font: inherit; }
    button {
      min-height: 36px;
      border: 1px solid var(--line);
      border-radius: 7px;
      background: #fff;
      color: var(--ink);
      padding: 0 12px;
      cursor: pointer;
      transition: border-color 140ms ease, box-shadow 140ms ease, background 140ms ease;
    }
    button:hover {
      border-color: #9fb3bb;
      box-shadow: 0 0 0 2px rgba(23, 107, 135, 0.08);
    }
    button.primary { border-color: var(--accent); background: var(--accent); color: #fff; font-weight: 800; }
    button.danger { border-color: rgba(184, 64, 58, 0.35); color: var(--danger); }
    input, textarea {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 7px;
      background: #fff;
      color: var(--ink);
      outline: none;
      transition: border-color 140ms ease, box-shadow 140ms ease;
    }
    input:focus, textarea:focus {
      border-color: var(--accent);
      box-shadow: 0 0 0 3px rgba(23, 107, 135, 0.13);
    }
    textarea {
      resize: vertical;
      min-height: 130px;
      padding: 10px;
      line-height: 1.45;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 13px;
    }
    h1, h2, h3, p { margin: 0; }
    h1 { margin-top: 3px; font-size: 22px; line-height: 1.1; }
    h2 { font-size: 18px; }
    h3 { font-size: 15px; }
    .app-shell {
      display: grid;
      grid-template-columns: 270px minmax(900px, 1fr);
      height: 100vh;
      overflow: hidden;
    }
    .sample-rail {
      display: flex;
      flex-direction: column;
      gap: 14px;
      min-height: 0;
      padding: 18px 14px;
      border-right: 1px solid var(--line);
      background: var(--rail);
    }
    .rail-head {
      display: flex;
      align-items: start;
      justify-content: space-between;
      gap: 12px;
    }
    .eyebrow {
      color: var(--accent);
      font-size: 11px;
      font-weight: 900;
      text-transform: uppercase;
      letter-spacing: .04em;
    }
    .status-pill {
      white-space: nowrap;
      border: 1px solid rgba(47, 143, 107, 0.3);
      border-radius: 999px;
      background: #eef8f2;
      color: #286f50;
      padding: 5px 9px;
      font-size: 12px;
      font-weight: 800;
    }
    .status-pill.saving {
      border-color: rgba(217, 138, 43, 0.4);
      background: #fff5e8;
      color: #9b5f17;
    }
    .status-pill.error {
      border-color: rgba(184, 64, 58, 0.35);
      background: #fff0ef;
      color: var(--danger);
    }
    .progress-panel {
      border: 1px solid rgba(216, 213, 204, 0.9);
      border-radius: var(--radius);
      background: rgba(255, 253, 248, 0.72);
      padding: 12px;
    }
    .priority-guide {
      border: 1px solid rgba(217, 138, 43, 0.34);
      border-radius: var(--radius);
      background: #fff8eb;
      color: #5f4630;
      padding: 10px 11px;
      font-size: 12px;
      line-height: 1.45;
    }
    .priority-guide strong {
      display: block;
      margin-bottom: 4px;
      color: #8a4f13;
      font-size: 13px;
    }
    .progress-top {
      display: flex;
      justify-content: space-between;
      color: var(--muted);
      font-size: 13px;
      font-weight: 800;
    }
    .progress-bar {
      height: 8px;
      margin-top: 10px;
      border-radius: 999px;
      background: #ded9cd;
      overflow: hidden;
    }
    #progressFill {
      width: 0%;
      height: 100%;
      border-radius: inherit;
      background: linear-gradient(90deg, var(--accent), var(--accent-2));
      transition: width 200ms ease;
    }
    .search-row input { height: 38px; padding: 0 11px; }
    .filter-row {
      display: grid;
      grid-template-columns: repeat(2, 1fr);
      gap: 8px;
    }
    .filter-row button {
      min-height: 34px;
      padding: 0 8px;
      color: var(--muted);
    }
    .filter-row button.active {
      border-color: var(--accent);
      background: #eaf4f7;
      color: var(--accent);
      font-weight: 900;
    }
    .sample-list {
      min-height: 0;
      overflow: auto;
      display: grid;
      grid-template-columns: 1fr;
      gap: 7px;
      padding-right: 2px;
    }
    .sample-tile {
      min-height: 76px;
      padding: 10px 11px;
      text-align: left;
      background: rgba(255, 253, 248, 0.78);
    }
    .sample-tile.priority {
      border-color: rgba(217, 138, 43, 0.42);
      background: #fffaf0;
    }
    .sample-tile.regular {
      background: rgba(255, 253, 248, 0.62);
    }
    .sample-tile.active {
      border-color: var(--accent);
      background: #f6ffff;
      box-shadow: 0 0 0 2px rgba(23, 107, 135, 0.16);
    }
    .tile-top {
      display: flex;
      align-items: baseline;
      justify-content: space-between;
      gap: 8px;
      min-width: 0;
    }
    .tile-id {
      display: block;
      min-width: 0;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      font-weight: 900;
      font-size: 16px;
      line-height: 1.1;
    }
    .tile-queue {
      flex: 0 0 auto;
      font-size: 12px;
      font-weight: 900;
    }
    .tile-queue.priority { color: var(--danger); }
    .tile-queue.regular { color: var(--muted); }
    .tile-meta {
      display: block;
      margin-top: 7px;
      color: var(--muted);
      font-size: 12px;
      line-height: 1.3;
    }
    .tile-status {
      display: block;
      margin-top: 5px;
      color: #53616b;
      font-size: 12px;
      font-weight: 800;
    }
    .pager-row {
      display: grid;
      grid-template-columns: 38px minmax(0, 1fr) 38px;
      gap: 7px;
      align-items: center;
    }
    .pager-row button {
      min-width: 38px;
      padding: 0;
      font-size: 18px;
      font-weight: 900;
    }
    .pager-row span {
      min-width: 0;
      overflow: hidden;
      text-align: center;
      color: var(--muted);
      font-size: 12px;
      font-weight: 800;
      white-space: nowrap;
      text-overflow: ellipsis;
    }
    .workspace {
      min-width: 0;
      min-height: 0;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .topbar {
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 16px;
      padding: 14px 16px;
      border-bottom: 1px solid var(--line);
      background: rgba(255, 253, 248, 0.82);
    }
    .sample-title {
      display: flex;
      align-items: center;
      gap: 12px;
      min-width: 0;
    }
    .sample-id {
      font-size: 28px;
      line-height: 1;
      font-weight: 950;
    }
    .sample-meta {
      margin-top: 4px;
      color: var(--muted);
      font-size: 13px;
      overflow: hidden;
      white-space: nowrap;
      text-overflow: ellipsis;
      max-width: 520px;
    }
    .icon-button {
      width: 38px;
      min-width: 38px;
      padding: 0;
      font-size: 26px;
      font-weight: 700;
    }
    .top-actions {
      display: flex;
      flex-wrap: wrap;
      justify-content: flex-end;
      gap: 8px;
    }
    .main-grid {
      min-height: 0;
      flex: 1;
      padding: 14px;
      display: grid;
      grid-template-columns: minmax(280px, 1fr) minmax(330px, 1.05fr) minmax(280px, .9fr);
      gap: 14px;
      overflow: hidden;
    }
    .pane {
      min-width: 0;
      min-height: 0;
      display: flex;
      flex-direction: column;
      border: 1px solid var(--line);
      border-radius: var(--radius);
      background: var(--panel);
      box-shadow: var(--shadow);
      overflow: hidden;
    }
    .pane-head {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 14px;
      border-bottom: 1px solid var(--line);
    }
    .pane-head p {
      margin-top: 4px;
      color: var(--muted);
      font-size: 12px;
    }
    .zoom-controls, .editor-actions {
      display: flex;
      align-items: center;
      gap: 7px;
      white-space: nowrap;
    }
    .zoom-controls button { min-width: 36px; padding: 0 9px; }
    .image-stage {
      flex: 1;
      min-height: 0;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: auto;
      background:
        linear-gradient(45deg, #f2f0ea 25%, transparent 25%),
        linear-gradient(-45deg, #f2f0ea 25%, transparent 25%),
        linear-gradient(45deg, transparent 75%, #f2f0ea 75%),
        linear-gradient(-45deg, transparent 75%, #f2f0ea 75%);
      background-size: 20px 20px;
      background-position: 0 0, 0 10px, 10px -10px, -10px 0px;
      padding: 18px;
    }
    .image-stage img {
      display: block;
      max-width: none;
      max-height: none;
      transform-origin: center center;
      user-select: none;
      -webkit-user-drag: none;
    }
    .editor-pane-body {
      min-height: 0;
      overflow: auto;
      padding: 12px;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .esmiles-card {
      border: 1px solid var(--line);
      border-radius: var(--radius);
      background: #fff;
      padding: 10px;
    }
    .esmiles-card-head {
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: center;
      margin-bottom: 8px;
    }
    .esmiles-card-head span {
      color: var(--muted);
      font-size: 12px;
      font-weight: 800;
    }
    .candidate-box, .reason-box {
      border: 1px solid #e4dfd3;
      border-radius: var(--radius);
      background: #fbfaf6;
      padding: 10px;
      color: var(--muted);
      line-height: 1.45;
      overflow-wrap: anywhere;
      font-size: 13px;
    }
    .candidate-box code {
      display: block;
      margin-top: 7px;
      padding: 8px;
      border-radius: 6px;
      background: #15202b;
      color: #e8f0f7;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 12px;
      white-space: pre-wrap;
    }
    .notes-box label {
      display: block;
      margin-bottom: 7px;
      color: var(--muted);
      font-size: 12px;
      font-weight: 800;
    }
    .render-stage {
      flex: 1;
      min-height: 260px;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: auto;
      padding: 14px;
      background: #fff;
    }
    .render-stage img {
      max-width: 100%;
      height: auto;
    }
    .inspect-body {
      min-height: 0;
      overflow: auto;
      padding: 12px;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .kv {
      display: grid;
      grid-template-columns: 88px minmax(0, 1fr);
      gap: 7px 10px;
      color: var(--muted);
      font-size: 13px;
      line-height: 1.35;
    }
    .kv strong { color: var(--ink); }
    .file-panel {
      border: 1px solid var(--line);
      border-radius: var(--radius);
      background: #fff;
      padding: 10px;
      color: var(--muted);
      font-size: 12px;
      line-height: 1.45;
      overflow-wrap: anywhere;
    }
    .export-panel {
      display: grid;
      gap: 9px;
      margin-top: auto;
    }
    .toast-stack {
      position: fixed;
      right: 18px;
      bottom: 18px;
      display: grid;
      gap: 8px;
      z-index: 9;
    }
    .toast {
      max-width: 360px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
      box-shadow: var(--shadow);
      padding: 10px 12px;
      font-size: 13px;
    }
    @media (max-width: 1100px) {
      body { min-width: 980px; }
      .app-shell { grid-template-columns: 270px minmax(710px, 1fr); }
      .main-grid {
        grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
        grid-template-rows: minmax(360px, .8fr) minmax(420px, 1fr);
        overflow: auto;
      }
      .image-pane, .editor-pane { min-height: 360px; }
      .inspect-pane { grid-column: 1 / -1; min-height: 420px; }
    }
  </style>
</head>
<body>
  <div class="app-shell">
    <aside class="sample-rail">
      <div class="rail-head">
        <div>
          <div class="eyebrow">Locus</div>
          <h1>E-SMILES 复查工作台</h1>
        </div>
        <div class="status-pill" id="saveStatus">启动中</div>
      </div>

      <div class="progress-panel">
        <div class="progress-top">
          <span id="progressText">0 / 4000 完成</span>
          <span id="changedText">0 修改 · 0 存疑</span>
        </div>
        <div class="progress-bar"><div id="progressFill"></div></div>
      </div>

      <div class="priority-guide">
        <strong>重点复查说明</strong>
        重点队列按风险分排序，建议优先人工复查；其他样本也在全量工作稿里，导出 zip 时会一起保留。
      </div>

      <div class="search-row">
        <input id="sampleSearch" type="search" placeholder="搜索文件名 / rank / 原因" autocomplete="off">
      </div>

      <div class="filter-row" id="filterRow">
        <button data-filter="all" class="active">全部</button>
        <button data-filter="priority">重点300</button>
        <button data-filter="todo">未完成</button>
        <button data-filter="changed">已修改</button>
        <button data-filter="flagged">存疑</button>
      </div>

      <div class="pager-row">
        <button id="listPrevBtn" title="上一页">‹</button>
        <span id="listPageInfo">载入中</span>
        <button id="listNextBtn" title="下一页">›</button>
      </div>

      <div class="sample-list" id="sampleList"></div>
    </aside>

    <main class="workspace">
      <header class="topbar">
        <div class="sample-title">
          <button class="icon-button" id="prevBtn" title="上一张">‹</button>
          <div>
            <div class="sample-id" id="sampleId">mol_0000.png</div>
            <div class="sample-meta" id="sampleMeta">载入中</div>
          </div>
          <button class="icon-button" id="nextBtn" title="下一张">›</button>
        </div>
        <div class="top-actions">
          <button id="restoreBtn">恢复初稿</button>
          <button id="useCandidateBtn">使用候选</button>
          <button id="saveBtn" class="primary">保存工作稿</button>
          <button id="markReviewedBtn" class="primary">完成这张</button>
          <button id="flagBtn">存疑</button>
        </div>
      </header>

      <section class="main-grid">
        <section class="pane image-pane">
          <div class="pane-head">
            <div>
              <h2>原图</h2>
              <p>左侧为比赛图片；滚轮缩放，按钮快速调整</p>
            </div>
            <div class="zoom-controls">
              <button id="zoomOutBtn">-</button>
              <span id="zoomLabel">100%</span>
              <button id="zoomInBtn">+</button>
              <button id="fitBtn">适配</button>
            </div>
          </div>
          <div class="image-stage" id="imageStage">
            <img id="originalImage" alt="original molecule" draggable="false">
          </div>
        </section>

        <section class="pane editor-pane">
          <div class="pane-head">
            <div>
              <h2>E-SMILES 编辑</h2>
              <p id="editState">改完会自动重画，Cmd+S 保存</p>
            </div>
            <div class="editor-actions">
              <button id="clearBtn">清空</button>
              <button id="copyBtn">复制</button>
            </div>
          </div>
          <div class="editor-pane-body">
            <article class="esmiles-card">
              <div class="esmiles-card-head">
                <h3>当前提交值</h3>
                <span id="changedBadge">未修改</span>
              </div>
              <textarea id="esmilesEditor" spellcheck="false"></textarea>
            </article>
            <div class="candidate-box">
              <strong>高分不同候选</strong>
              <div id="candidateMeta">无候选</div>
              <code id="candidateText"></code>
            </div>
            <div class="reason-box" id="reasonBox"></div>
            <div class="notes-box">
              <label for="sampleNote">这张图的备注</label>
              <textarea id="sampleNote" rows="3" placeholder="例如：环系疑似错；手性需再查"></textarea>
            </div>
          </div>
        </section>

        <section class="pane inspect-pane">
          <div class="pane-head">
            <div>
              <h2>RDKit 实时核查</h2>
              <p id="rdkitState">等待输入</p>
            </div>
            <button id="buildZipBtn">导出工作 zip</button>
          </div>
          <div class="render-stage">
            <img id="renderedImage" alt="rendered molecule">
          </div>
          <div class="inspect-body">
            <div class="kv">
              <strong>队列</strong><span id="scopeValue">-</span>
              <strong>风险分</strong><span id="riskValue">-</span>
              <strong>排名</strong><span id="rankValue">-</span>
              <strong>页位</strong><span id="pageValue">-</span>
              <strong>原分</strong><span id="scoreValue">-</span>
              <strong>最佳候选</strong><span id="bestValue">-</span>
            </div>
            <div class="file-panel">
              <strong>保存文件</strong><br>
              working_submission.csv：全量工作稿<br>
              review_overrides.csv：只包含你改过的行<br>
              working_submission_package.zip：点击导出后生成
            </div>
            <div class="export-panel">
              <button id="downloadCsvBtn">下载 CSV</button>
              <button id="downloadZipBtn">下载 zip</button>
              <button id="showChangedBtn">只看已修改</button>
              <button id="showAllBtn">回到全部</button>
            </div>
          </div>
        </section>
      </section>
    </main>
  </div>

  <div class="toast-stack" id="toastStack"></div>

  <script>
    const state = {
      items: [],
      filtered: [],
      activeIndex: 0,
      active: null,
      filter: "all",
      query: "",
      listPage: 0,
      pageSize: 80,
      review: {},
      zoom: 1,
      renderTimer: null,
      autoSaveTimer: null,
      reviewSaveTimer: null,
      dirty: false,
    };
    const els = {};
    const reviewKey = "locus-review-editor-v2";

    document.addEventListener("DOMContentLoaded", () => {
      init().catch((error) => {
        setSaveStatus("启动失败", "error");
        toast(`启动失败：${error}`);
      });
    });

    async function init() {
      bindElements();
      bindEvents();
      await loadReviewState();
      await loadItems();
    }

    function bindElements() {
      [
        "saveStatus", "progressText", "changedText", "progressFill", "sampleSearch", "filterRow", "sampleList",
        "listPrevBtn", "listPageInfo", "listNextBtn",
        "prevBtn", "nextBtn", "sampleId", "sampleMeta", "restoreBtn", "useCandidateBtn", "saveBtn",
        "markReviewedBtn", "flagBtn", "zoomOutBtn", "zoomInBtn", "fitBtn", "zoomLabel", "imageStage",
        "originalImage", "clearBtn", "copyBtn", "editState", "changedBadge", "esmilesEditor",
        "candidateMeta", "candidateText", "reasonBox", "sampleNote", "rdkitState", "buildZipBtn",
        "renderedImage", "scopeValue", "riskValue", "rankValue", "pageValue", "scoreValue", "bestValue",
        "downloadCsvBtn", "downloadZipBtn", "showChangedBtn", "showAllBtn", "toastStack"
      ].forEach((id) => els[id] = document.getElementById(id));
    }

    function bindEvents() {
      els.sampleSearch.addEventListener("input", () => {
        state.query = els.sampleSearch.value.trim();
        state.listPage = 0;
        renderSampleList();
      });
      els.filterRow.addEventListener("click", (event) => {
        const button = event.target.closest("button[data-filter]");
        if (!button) return;
        state.filter = button.dataset.filter;
        state.listPage = 0;
        renderFilters();
        renderSampleList();
        els.sampleList.scrollTop = 0;
      });
      els.listPrevBtn.addEventListener("click", () => changeListPage(-1));
      els.listNextBtn.addEventListener("click", () => changeListPage(1));
      els.prevBtn.addEventListener("click", () => moveSample(-1));
      els.nextBtn.addEventListener("click", () => moveSample(1));
      els.restoreBtn.addEventListener("click", restoreOriginal);
      els.useCandidateBtn.addEventListener("click", useCandidate);
      els.saveBtn.addEventListener("click", saveActive);
      els.markReviewedBtn.addEventListener("click", toggleReviewed);
      els.flagBtn.addEventListener("click", toggleFlagged);
      els.clearBtn.addEventListener("click", () => setEditorValue(""));
      els.copyBtn.addEventListener("click", copyEsmiles);
      els.sampleNote.addEventListener("input", () => {
        reviewForActive().note = els.sampleNote.value;
        saveReviewState();
        renderSampleList();
      });
      els.esmilesEditor.addEventListener("input", () => {
        state.dirty = true;
        setSaveStatus("等待自动保存", "saving");
        updateChangedBadge();
        scheduleRender();
        scheduleAutoSave();
      });
      els.zoomOutBtn.addEventListener("click", () => setZoom(state.zoom / 1.18));
      els.zoomInBtn.addEventListener("click", () => setZoom(state.zoom * 1.18));
      els.fitBtn.addEventListener("click", fitImage);
      els.originalImage.addEventListener("load", () => {
        fitImage();
        setTimeout(fitImage, 80);
        setTimeout(fitImage, 280);
      });
      els.imageStage.addEventListener("wheel", onWheel, {passive: false});
      els.buildZipBtn.addEventListener("click", buildZip);
      els.downloadCsvBtn.addEventListener("click", downloadCsv);
      els.downloadZipBtn.addEventListener("click", downloadZip);
      els.showChangedBtn.addEventListener("click", () => { state.filter = "changed"; state.listPage = 0; renderFilters(); renderSampleList(); els.sampleList.scrollTop = 0; });
      els.showAllBtn.addEventListener("click", () => { state.filter = "all"; state.listPage = 0; state.query = ""; els.sampleSearch.value = ""; renderFilters(); renderSampleList(); els.sampleList.scrollTop = 0; });
      document.addEventListener("keydown", (event) => {
        if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === "s") {
          event.preventDefault();
          saveActive();
          return;
        }
        if (isTyping(event.target)) return;
        if (event.key === "ArrowLeft") moveSample(-1);
        if (event.key === "ArrowRight") moveSample(1);
        if (event.key.toLowerCase() === "d") toggleReviewed();
        if (event.key.toLowerCase() === "f") toggleFlagged();
      });
    }

    async function loadItems() {
      setSaveStatus("载入中", "saving");
      const response = await fetch("/api/items");
      state.items = await response.json();
      state.filtered = state.items.slice();
      setSaveStatus("已自动保存", "");
      renderFilters();
      renderSampleList();
      renderProgress();
      setActiveByIndex(0);
    }

    async function loadReviewState() {
      let localState = {};
      try {
        localState = JSON.parse(localStorage.getItem(reviewKey) || "{}");
      } catch {
        localState = {};
      }
      try {
        const response = await fetch("/api/review_state");
        const data = await response.json();
        state.review = data.ok && data.review && Object.keys(data.review).length ? data.review : localState;
      } catch {
        state.review = localState;
      }
    }

    function saveReviewState() {
      localStorage.setItem(reviewKey, JSON.stringify(state.review));
      renderProgress();
      clearTimeout(state.reviewSaveTimer);
      state.reviewSaveTimer = setTimeout(saveReviewStateToDisk, 500);
    }

    async function saveReviewStateToDisk() {
      try {
        await fetch("/api/review_state", {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify({review: state.review})
        });
      } catch {
        setSaveStatus("备注保存失败", "error");
      }
    }

    function reviewFor(id) {
      if (!state.review[id]) state.review[id] = {reviewed: false, flagged: false, note: ""};
      return state.review[id];
    }

    function reviewForActive() {
      return reviewFor(state.active?.file_name || "");
    }

    function renderFilters() {
      els.filterRow.querySelectorAll("button").forEach((button) => {
        button.classList.toggle("active", button.dataset.filter === state.filter);
      });
    }

    function filteredItems() {
      const q = state.query.toLowerCase();
      return state.items.filter((item) => {
        const review = reviewFor(item.file_name);
        const changed = item.e_smiles !== item.original_e_smiles;
        const priority = item.is_priority_review === "true";
        if (state.filter === "priority" && !priority) return false;
        if (state.filter === "todo" && review.reviewed) return false;
        if (state.filter === "changed" && !changed) return false;
        if (state.filter === "flagged" && !review.flagged) return false;
        if (!q) return true;
        return [
          item.file_name, item.rank, item.page, item.slot, item.risk, item.reasons, item.reasons_zh,
          item.priority_label, item.priority_note,
          item.e_smiles, item.diff_candidate
        ].join(" ").toLowerCase().includes(q);
      });
    }

    function renderSampleList() {
      state.filtered = filteredItems();
      const total = state.filtered.length;
      const totalPages = Math.max(1, Math.ceil(total / state.pageSize));
      state.listPage = Math.max(0, Math.min(state.listPage, totalPages - 1));
      const start = state.listPage * state.pageSize;
      const end = Math.min(total, start + state.pageSize);
      const visibleItems = state.filtered.slice(start, end);
      const fragment = document.createDocumentFragment();
      if (!visibleItems.length) {
        const empty = document.createElement("div");
        empty.className = "file-panel";
        empty.textContent = "没有符合条件的样本。可以切回“全部”或换个搜索词。";
        fragment.appendChild(empty);
      }
      for (const item of visibleItems) {
        const review = reviewFor(item.file_name);
        const changed = item.e_smiles !== item.original_e_smiles;
        const priority = item.is_priority_review === "true";
        const button = document.createElement("button");
        button.className = "sample-tile";
        button.classList.toggle("priority", priority);
        button.classList.toggle("regular", !priority);
        button.classList.toggle("active", state.active?.file_name === item.file_name);
        const queueText = priority ? `重点 #${item.rank || "-"}` : "全量";
        const statusText = [
          changed ? "已修改" : "未修改",
          review.reviewed ? "已完成" : "未完成",
          review.flagged ? "存疑" : ""
        ].filter(Boolean).join(" · ");
        const meta = priority
          ? `风险 ${fmt(item.risk)} · 第 ${item.rank || "-"} / ${item.risk_total || "300"} · p${item.page || "-"}-${item.slot || "-"}`
          : `未进重点300 · 全量可搜索/可修改`;
        button.innerHTML = `
          <span class="tile-top">
            <span class="tile-id">${escapeHtml(item.file_name.replace(".png", ""))}</span>
            <span class="tile-queue ${priority ? "priority" : "regular"}">${escapeHtml(queueText)}</span>
          </span>
          <span class="tile-meta">${escapeHtml(meta)}</span>
          <span class="tile-status">${escapeHtml(statusText)}</span>
        `;
        button.addEventListener("click", () => navigateToFile(item.file_name));
        fragment.appendChild(button);
      }
      els.sampleList.replaceChildren(fragment);
      els.listPageInfo.textContent = total
        ? `${start + 1}-${end} / ${total}`
        : "0 / 0";
      els.listPrevBtn.disabled = state.listPage <= 0;
      els.listNextBtn.disabled = state.listPage >= totalPages - 1;
      renderProgress();
    }

    function changeListPage(delta) {
      const totalPages = Math.max(1, Math.ceil(state.filtered.length / state.pageSize));
      state.listPage = Math.max(0, Math.min(state.listPage + delta, totalPages - 1));
      renderSampleList();
      els.sampleList.scrollTop = 0;
    }

    function renderProgress() {
      const total = state.items.length || 1;
      const priorityItems = state.items.filter((item) => item.is_priority_review === "true");
      const priorityTotal = priorityItems.length || 1;
      const reviewed = state.items.filter((item) => reviewFor(item.file_name).reviewed).length;
      const reviewedPriority = priorityItems.filter((item) => reviewFor(item.file_name).reviewed).length;
      const changed = state.items.filter((item) => item.e_smiles !== item.original_e_smiles).length;
      const flagged = state.items.filter((item) => reviewFor(item.file_name).flagged).length;
      els.progressText.textContent = `重点 ${reviewedPriority}/${priorityTotal} · 全量 ${reviewed}/${total}`;
      els.changedText.textContent = `${changed} 修改 · ${flagged} 存疑`;
      els.progressFill.style.width = `${Math.round((reviewedPriority / priorityTotal) * 100)}%`;
    }

    function setActive(fileName) {
      const index = state.filtered.findIndex((item) => item.file_name === fileName);
      if (index >= 0) setActiveByIndex(index);
    }

    async function navigateToFile(fileName) {
      if (state.dirty) {
        const ok = await saveActive({quiet: true});
        if (!ok) return;
      }
      setActive(fileName);
    }

    function setActiveByIndex(index) {
      if (!state.filtered.length) return;
      state.activeIndex = Math.max(0, Math.min(index, state.filtered.length - 1));
      state.listPage = Math.floor(state.activeIndex / state.pageSize);
      state.active = state.filtered[state.activeIndex];
      state.dirty = false;
      clearTimeout(state.autoSaveTimer);
      const item = state.active;
      const review = reviewFor(item.file_name);
      const priority = item.is_priority_review === "true";
      els.sampleId.textContent = item.file_name;
      els.sampleMeta.textContent = priority
        ? `重点复查 · 第 ${item.rank || "-"} / ${item.risk_total || "300"} · page ${item.page || "-"}-${item.slot || "-"} · 风险分 ${fmt(item.risk)} · 原分 ${item.orig_score || "-"}`
        : `全量样本 · 未进入前 ${item.risk_total || "300"} 重点风险队列 · 导出仍包含这张`;
      els.originalImage.src = `/preview/${encodeURIComponent(item.file_name)}`;
      els.esmilesEditor.value = item.e_smiles || "";
      els.sampleNote.value = review.note || "";
      els.candidateMeta.textContent = item.diff_candidate ? `${item.diff_source_zh || item.diff_source || "候选"} · 分数 ${item.diff_score || "-"}` : "无候选";
      els.candidateText.textContent = item.diff_candidate || "";
      els.reasonBox.innerHTML = `<strong>${priority ? "中文风险原因" : "队列说明"}</strong><br>${lineBreaks(item.reasons_zh || "无")}`;
      els.scopeValue.textContent = priority ? `重点复查前 ${item.risk_total || "300"}` : "全量样本";
      els.riskValue.textContent = priority ? fmt(item.risk) : "未列入";
      els.rankValue.textContent = priority ? `${item.rank || "-"} / ${item.risk_total || "300"}` : "未进前300";
      els.pageValue.textContent = priority ? `${item.page || "-"}-${item.slot || "-"}` : "-";
      els.scoreValue.textContent = item.orig_score || "-";
      els.bestValue.textContent = `${item.best_source_zh || item.best_source || "-"} ${item.best_score || ""}`;
      els.markReviewedBtn.textContent = review.reviewed ? "取消完成" : "完成这张";
      els.flagBtn.textContent = review.flagged ? "取消存疑" : "存疑";
      updateChangedBadge();
      renderSampleList();
      requestAnimationFrame(fitImage);
      setTimeout(fitImage, 120);
      scheduleRender(0);
    }

    function moveSample(delta) {
      navigateByIndex(state.activeIndex + delta);
    }

    async function navigateByIndex(index) {
      if (state.dirty) {
        const ok = await saveActive({quiet: true});
        if (!ok) return;
      }
      setActiveByIndex(index);
    }

    function setEditorValue(value) {
      els.esmilesEditor.value = value;
      state.dirty = true;
      setSaveStatus("等待自动保存", "saving");
      updateChangedBadge();
      scheduleRender(0);
      scheduleAutoSave();
    }

    function restoreOriginal() {
      if (!state.active) return;
      setEditorValue(state.active.original_e_smiles || "");
    }

    function useCandidate() {
      if (!state.active?.diff_candidate) {
        toast("没有记录不同候选");
        return;
      }
      setEditorValue(state.active.diff_candidate);
    }

    async function copyEsmiles() {
      await navigator.clipboard.writeText(els.esmilesEditor.value);
      toast("已复制 E-SMILES");
    }

    function toggleReviewed() {
      if (!state.active) return;
      const review = reviewForActive();
      review.reviewed = !review.reviewed;
      saveReviewState();
      setActiveByIndex(state.activeIndex);
    }

    function toggleFlagged() {
      if (!state.active) return;
      const review = reviewForActive();
      review.flagged = !review.flagged;
      saveReviewState();
      setActiveByIndex(state.activeIndex);
    }

    function updateChangedBadge() {
      if (!state.active) return;
      const value = els.esmilesEditor.value.trim();
      const changed = value !== (state.active.original_e_smiles || "");
      els.changedBadge.textContent = state.dirty ? "等待自动保存" : (changed ? "已修改" : "未修改");
      els.changedBadge.style.color = state.dirty ? "var(--warm)" : (changed ? "var(--accent-2)" : "var(--muted)");
      els.editState.textContent = state.dirty ? "改动会自动保存，也可 Cmd+S 立即保存" : "改完会自动重画并自动保存";
    }

    function scheduleRender(delay = 700) {
      clearTimeout(state.renderTimer);
      state.renderTimer = setTimeout(renderCurrent, delay);
    }

    function scheduleAutoSave(delay = 900) {
      clearTimeout(state.autoSaveTimer);
      state.autoSaveTimer = setTimeout(() => saveActive({quiet: true}), delay);
    }

    async function renderCurrent() {
      const esmi = els.esmilesEditor.value;
      els.rdkitState.textContent = "渲染中";
      const response = await fetch("/api/render", {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify({e_smiles: esmi})
      });
      const blob = await response.blob();
      const old = els.renderedImage.src;
      els.renderedImage.src = URL.createObjectURL(blob);
      if (old.startsWith("blob:")) URL.revokeObjectURL(old);
      els.rdkitState.textContent = "已根据当前 E-SMILES 重画";
    }

    async function saveActive(options = {}) {
      const quiet = Boolean(options.quiet);
      if (!state.active) return true;
      clearTimeout(state.autoSaveTimer);
      setSaveStatus("保存中", "saving");
      const fileName = state.active.file_name;
      const value = els.esmilesEditor.value.trim();
      const response = await fetch("/api/save", {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify({file_name: fileName, e_smiles: value})
      });
      const data = await response.json();
      if (!data.ok) {
        setSaveStatus("保存失败", "error");
        toast(data.error || "保存失败");
        return false;
      }
      const global = state.items.find((item) => item.file_name === fileName);
      if (global) {
        global.e_smiles = data.e_smiles;
        if (state.active?.file_name === fileName) Object.assign(state.active, global);
      }
      if (state.active?.file_name === fileName) state.dirty = false;
      setSaveStatus("已自动保存", "");
      updateChangedBadge();
      renderSampleList();
      if (!quiet) toast("已保存到 CSV");
      return true;
    }

    async function buildZip() {
      if (state.dirty) {
        const ok = await saveActive({quiet: true});
        if (!ok) return;
      }
      await saveReviewStateToDisk();
      setSaveStatus("导出中", "saving");
      const response = await fetch("/api/build_zip", {method: "POST"});
      const data = await response.json();
      if (data.ok) {
        setSaveStatus("已自动保存", "");
        toast(`已生成 ${data.path}`);
      } else {
        setSaveStatus("导出失败", "error");
        toast(data.error || "导出失败");
      }
    }

    async function downloadCsv() {
      if (state.dirty) {
        const ok = await saveActive({quiet: true});
        if (!ok) return;
      }
      await saveReviewStateToDisk();
      window.location.href = "/download/submission.csv";
    }

    async function downloadZip() {
      if (state.dirty) {
        const ok = await saveActive({quiet: true});
        if (!ok) return;
      }
      await saveReviewStateToDisk();
      setSaveStatus("导出中", "saving");
      const response = await fetch("/api/build_zip", {method: "POST"});
      const data = await response.json();
      if (!data.ok) {
        setSaveStatus("导出失败", "error");
        toast(data.error || "导出失败");
        return;
      }
      setSaveStatus("已自动保存", "");
      window.location.href = "/download/working_submission_package.zip";
    }

    function setSaveStatus(text, mode) {
      els.saveStatus.textContent = text;
      els.saveStatus.classList.toggle("saving", mode === "saving");
      els.saveStatus.classList.toggle("error", mode === "error");
    }

    function setZoom(value) {
      state.zoom = Math.max(0.04, Math.min(5, value));
      if (els.originalImage.naturalWidth) {
        els.originalImage.style.width = `${Math.round(els.originalImage.naturalWidth * state.zoom)}px`;
        els.originalImage.style.height = "auto";
        els.originalImage.style.transform = "none";
      } else {
        els.originalImage.style.transform = `scale(${state.zoom})`;
      }
      els.zoomLabel.textContent = `${Math.round(state.zoom * 100)}%`;
    }

    function fitImage() {
      const img = els.originalImage;
      const stage = els.imageStage;
      if (!img.naturalWidth || !stage.clientWidth) {
        setZoom(1);
        return;
      }
      const pad = 36;
      const rect = stage.getBoundingClientRect();
      const availableWidth = Math.max(120, Math.min(stage.clientWidth, rect.width) - pad);
      const availableHeight = Math.max(120, Math.min(stage.clientHeight, rect.height) - pad);
      const scale = Math.min(
        availableWidth / img.naturalWidth,
        availableHeight / img.naturalHeight,
        1
      );
      setZoom(scale);
      stage.scrollTo({top: 0, left: 0, behavior: "auto"});
    }

    function onWheel(event) {
      event.preventDefault();
      setZoom(state.zoom * (event.deltaY < 0 ? 1.12 : 1 / 1.12));
    }

    function isTyping(target) {
      const tag = target?.tagName?.toLowerCase();
      return tag === "input" || tag === "textarea";
    }

    function escapeHtml(value) {
      return String(value || "").replace(/[&<>"']/g, (c) => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[c]));
    }

    function lineBreaks(value) {
      return escapeHtml(value).replace(/\n/g, "<br>");
    }

    function fmt(value) {
      const n = Number(value);
      return Number.isFinite(n) ? n.toFixed(1) : "-";
    }

    function toast(message) {
      const el = document.createElement("div");
      el.className = "toast";
      el.textContent = message;
      els.toastStack.appendChild(el);
      setTimeout(() => el.remove(), 2600);
    }
  </script>
</body>
</html>
"""


class Handler(BaseHTTPRequestHandler):
    state: ReviewState

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def send_bytes(self, body: bytes, content_type: str, status: int = 200) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_download(self, path: Path, download_name: str, content_type: str) -> None:
        body = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Disposition", f'attachment; filename="{download_name}"')
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_json(self, value: Any, status: int = 200) -> None:
        self.send_bytes(json.dumps(value, ensure_ascii=False).encode("utf-8"), "application/json; charset=utf-8", status)

    def read_json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length <= 0:
            return {}
        return json.loads(self.rfile.read(length).decode("utf-8"))

    def do_GET(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        if path == "/":
            self.send_bytes(WORKBENCH_HTML.encode("utf-8"), "text/html; charset=utf-8")
            return
        if path == "/api/items":
            with self.state.lock:
                self.send_json([self.state.item_payload(row) for row in self.state.risk_rows])
            return
        if path == "/api/review_state":
            with self.state.lock:
                self.send_json({"ok": True, "review": self.state.read_review_state()})
            return
        if path == "/download/submission.csv":
            with self.state.lock:
                self.send_download(self.state.working_csv, "submission.csv", "text/csv; charset=utf-8")
            return
        if path == "/download/working_submission_package.zip":
            with self.state.lock:
                if not self.state.package_zip.exists():
                    self.state.build_zip()
                self.send_download(self.state.package_zip, "working_submission_package.zip", "application/zip")
            return
        if path.startswith("/preview/"):
            name = Path(urllib.parse.unquote(path.removeprefix("/preview/"))).name
            src = self.state.pic_dir / name
            if not src.exists() or src.parent != self.state.pic_dir:
                self.send_json({"ok": False, "error": "image not found"}, status=404)
                return
            cache = self.state.preview_dir / f"{src.stem}.png"
            try:
                body = image_preview_png(src, cache)
            except Exception:
                body = src.read_bytes()
                content_type = "image/png" if src.suffix.lower() == ".png" else "image/jpeg"
                self.send_bytes(body, content_type)
                return
            self.send_bytes(body, "image/png")
            return
        if path.startswith("/original/"):
            name = Path(urllib.parse.unquote(path.removeprefix("/original/"))).name
            src = self.state.pic_dir / name
            if not src.exists() or src.parent != self.state.pic_dir:
                self.send_json({"ok": False, "error": "image not found"}, status=404)
                return
            content_type = "image/png" if src.suffix.lower() == ".png" else "image/jpeg"
            self.send_bytes(src.read_bytes(), content_type)
            return
        self.send_json({"ok": False, "error": "not found"}, status=404)

    def do_POST(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path == "/api/render":
            data = self.read_json()
            png = render_esmi_png(str(data.get("e_smiles", "")))
            self.send_bytes(png, "image/png")
            return
        if parsed.path == "/api/save":
            data = self.read_json()
            name = str(data.get("file_name", ""))
            esmi = ensure_esmi(str(data.get("e_smiles", "")))
            with self.state.lock:
                if name not in self.state.values:
                    self.send_json({"ok": False, "error": "unknown file_name"}, status=400)
                    return
                self.state.values[name] = esmi
                self.state.save_files()
                self.send_json(
                    {
                        "ok": True,
                        "file_name": name,
                        "e_smiles": esmi,
                        "changed": esmi != self.state.original_values.get(name, ""),
                        "working_csv": str(self.state.working_csv),
                        "overrides_csv": str(self.state.overrides_csv),
                    }
                )
            return
        if parsed.path == "/api/review_state":
            data = self.read_json()
            review = data.get("review", {})
            if not isinstance(review, dict):
                self.send_json({"ok": False, "error": "review must be an object"}, status=400)
                return
            with self.state.lock:
                self.state.save_review_state(review)
                self.send_json({"ok": True, "path": str(self.state.review_state_json)})
            return
        if parsed.path == "/api/build_zip":
            with self.state.lock:
                try:
                    self.state.build_zip()
                except Exception as exc:
                    self.send_json({"ok": False, "error": repr(exc)}, status=500)
                    return
                self.send_json({"ok": True, "path": str(self.state.package_zip)})
            return
        self.send_json({"ok": False, "error": "not found"}, status=404)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pic-dir", required=True, type=Path)
    parser.add_argument("--submission-csv", required=True, type=Path)
    parser.add_argument("--risk-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--meta-path", type=Path)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()

    state = ReviewState(args)
    Handler.state = state
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(
        json.dumps(
            {
                "event": "ready",
                "url": f"http://{args.host}:{args.port}/",
                "items": len(state.risk_rows),
                "working_csv": str(state.working_csv),
                "overrides_csv": str(state.overrides_csv),
            },
            ensure_ascii=False,
        ),
        flush=True,
    )
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
