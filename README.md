# Locus

Locus is an open-source workflow for optical chemical structure recognition (OCSR) and E-SMILES submission packaging. It was prepared for [Bohrium Science Data Annotation Competition, Track 1](https://www.bohrium.com/competitions/53859761357?tab=introduce), where each molecule image must be converted into an Extended-SMILES result.

The project is designed as a practical competition toolkit rather than a single neural network. It batches recognition requests, generates safer image variants, merges multiple model outputs, renders predictions back into molecule images, and builds the final `submission.csv` + `meta.md` zip expected by the competition platform.

This repository does **not** include the official competition images, private submission files, generated zip packages, third-party model weights, or proprietary MolParser service code.

## What It Does

- Sends molecule images to the public MolParser OCSR endpoint and stores resumable CSV results.
- Creates preprocessing variants such as crop, padding, upscaling, binary cleanup, and color-foreground extraction for difficult images.
- Optionally runs local MolScribe or DECIMER helpers when they are installed separately.
- Normalizes model outputs into the competition's E-SMILES style.
- Uses RDKit to check whether the base SMILES is parseable.
- Renders predicted structures back to images and produces visual review sheets for human audit.
- Builds a platform-ready zip containing only `submission.csv` and `meta.md`.

## Repository Layout

```text
Locus/
├── tools/
│   ├── run_molparser_api.py          # Batch MolParser calls on original images
│   ├── run_molparser_variants.py     # Batch MolParser calls on preprocessed image variants
│   ├── run_molscribe.py              # Optional MolScribe fallback runner
│   ├── run_decimer.py                # Optional DECIMER fallback runner
│   ├── audit_v1_render_compare.py    # RDKit rendering and visual-risk audit
│   └── make_submission.py            # Merge candidates and build the final zip
├── examples/
│   └── overrides.example.csv         # Public example of the manual correction format
├── requirements.txt
├── LICENSE
└── README.md
```

## Requirements

Core workflow:

- Python 3.10 or newer is recommended.
- Pillow
- NumPy
- OpenCV
- scikit-image
- RDKit

Optional helpers:

- MolScribe, installed from its upstream project and used with a separately downloaded checkpoint.
- DECIMER, installed from its upstream project.
- A working internet connection for MolParser public OCSR calls.

Please review and follow the terms of use and licenses of any third-party model or web service you call. This repository only contains orchestration code.

## Installation

Clone the repository:

```bash
git clone https://github.com/Shawn-TV/Locus.git
cd Locus
```

Create a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If RDKit is difficult to install with pip on your platform, a Conda environment is often easier:

```bash
conda create -n locus python=3.11 rdkit pillow numpy opencv scikit-image -c conda-forge
conda activate locus
```

Install optional tools only if you plan to use them. For MolScribe and DECIMER, follow their upstream installation instructions and keep model checkpoints outside this repository.

## Data Preparation

Place the competition files in a local data directory. A typical local layout is:

```text
data/
├── pic/
│   ├── mol_0001.png
│   ├── mol_0002.png
│   └── ...
└── submission_template/
    └── submission.csv
```

The `data/` directory is ignored by Git on purpose, because the competition images and templates are governed by the competition organizer's rules.

Create a working directory for generated CSVs:

```bash
mkdir -p work build
```

## Quick Start

Run MolParser on the original images:

```bash
python tools/run_molparser_api.py \
  --input-dir data/pic \
  --output work/molparser_original.csv \
  --workers 6 \
  --resume
```

Run MolParser again on safer image variants, especially for low-confidence rows:

```bash
python tools/run_molparser_variants.py \
  --input-dir data/pic \
  --from-csv work/molparser_original.csv \
  --threshold 0.95 \
  --include-aggressive \
  --output work/molparser_variants.csv
```

If you are not using MolScribe, create an empty fallback CSV:

```bash
printf 'file_name,e_smiles,raw_smiles,confidence,error\n' > work/molscribe.csv
```

If you have MolScribe installed and a checkpoint available, run it as an optional fallback:

```bash
python tools/run_molscribe.py \
  --ckpt /path/to/molscribe.ckpt \
  --input-dir data/pic \
  --output work/molscribe.csv \
  --device cpu \
  --resume
```

Build an initial submission zip:

```bash
python tools/make_submission.py \
  --pic-dir data/pic \
  --template-csv data/submission_template/submission.csv \
  --molparser-csv work/molparser_original.csv work/molparser_variants.csv \
  --molscribe-csv work/molscribe.csv \
  --version V1 \
  --code-repo https://github.com/Shawn-TV/Locus \
  --out-dir build/V1 \
  --zip-path build/V1.zip
```

Audit the result by rendering the predicted E-SMILES back into images:

```bash
python tools/audit_v1_render_compare.py \
  --v1-zip build/V1.zip \
  --pic-dir data/pic \
  --out-dir reports/V1-audit \
  --molparser-csv work/molparser_original.csv work/molparser_variants.csv \
  --molscribe-csv work/molscribe.csv \
  --review-top 180
```

Review the generated sheets in `reports/V1-audit/review_sheets/`. If you need manual corrections, create a private override CSV using the schema in `examples/overrides.example.csv`, then rebuild:

```bash
python tools/make_submission.py \
  --pic-dir data/pic \
  --template-csv data/submission_template/submission.csv \
  --molparser-csv work/molparser_original.csv work/molparser_variants.csv \
  --molscribe-csv work/molscribe.csv \
  --overrides-csv /path/to/private_overrides.csv \
  --version V2 \
  --code-repo https://github.com/Shawn-TV/Locus \
  --out-dir build/V2 \
  --zip-path build/V2.zip
```

## Output

The final zip contains:

```text
submission.csv
meta.md
```

`submission.csv` uses the required columns:

```csv
file_name,e_smiles
mol_0001.png,CCO<sep>
```

`meta.md` records the model/tool usage, the code repository link, and a short method description.

## Notes on E-SMILES

Locus preserves MolParser E-SMILES captions whenever possible, because the competition expects Extended-SMILES syntax rather than plain SMILES only. The merger script also ensures populated rows include `<sep>` and keeps E-SMILES extension tags such as atom labels.

Some abbreviated groups or Markush labels can require human review. The audit script is built for that loop: render, compare, inspect, correct, and rebuild.

## Third-Party Tools and Licenses

Locus is released under the MIT License. Third-party tools remain under their own licenses and terms:

- MolParser public OCSR service: provided by DP Technology.
- MolScribe: optional open-source model/checkpoint from its upstream project.
- RDKit: open-source cheminformatics toolkit.
- DECIMER: optional OCSR helper from its upstream project.
- Bohrium competition data: provided by the competition organizer and not redistributed here.

## Limitations

This toolkit helps organize an OCSR competition workflow, but it does not guarantee chemical correctness. Low-confidence predictions, Markush structures, stereochemistry, bridged rings, abbreviations, and scanned labels still need careful human review.

---

# Locus 中文说明

Locus 是一个开源的分子结构图识别（OCSR）和 E-SMILES 提交包生成流程。它是为 [Bohrium 科学数据标注大赛赛道 1](https://www.bohrium.com/competitions/53859761357?tab=introduce) 准备的：赛题要求把每张分子图片识别成 Extended-SMILES，并按平台格式提交。

这个项目不是单一神经网络，而是一个实用的比赛工具包。它可以批量调用识别服务、生成更容易识别的图片预处理版本、合并多个模型结果、把预测结果重新画成分子图用于人工核验，最后生成比赛要求的 `submission.csv` + `meta.md` 压缩包。

本仓库**不包含**官方比赛图片、私人提交结果、生成好的 zip 包、第三方模型权重，也不包含 MolParser 的专有服务代码。

## 它能做什么

- 批量把分子图片发送给 MolParser public OCSR 接口，并保存可断点续跑的 CSV 结果。
- 为疑难图片生成预处理版本，例如裁剪、加白边、放大、黑白化、彩色前景提取等。
- 在你单独安装后，可选运行本地 MolScribe 或 DECIMER 作为补充识别来源。
- 把模型输出整理成比赛需要的 E-SMILES 风格。
- 使用 RDKit 检查底层 SMILES 是否可解析。
- 把预测结构重新渲染成分子图片，生成方便人工审阅的对比图。
- 打包生成平台可提交的 zip，里面只包含 `submission.csv` 和 `meta.md`。

## 仓库结构

```text
Locus/
├── tools/
│   ├── run_molparser_api.py          # 对原图批量调用 MolParser
│   ├── run_molparser_variants.py     # 对预处理图片版本批量调用 MolParser
│   ├── run_molscribe.py              # 可选的 MolScribe 补充识别脚本
│   ├── run_decimer.py                # 可选的 DECIMER 补充识别脚本
│   ├── audit_v1_render_compare.py    # RDKit 回画与视觉风险审计
│   └── make_submission.py            # 合并结果并生成最终提交包
├── examples/
│   └── overrides.example.csv         # 公开的人工修正表格式示例
├── requirements.txt
├── LICENSE
└── README.md
```

## 环境要求

核心流程需要：

- 推荐 Python 3.10 或更新版本。
- Pillow
- NumPy
- OpenCV
- scikit-image
- RDKit

可选工具：

- MolScribe：需要按上游项目说明安装，并自行下载 checkpoint。
- DECIMER：需要按上游项目说明单独安装。
- MolParser public OCSR 调用需要联网。

请自行遵守所有第三方模型和服务的使用条款与许可证。本仓库只提供流程编排代码。

## 安装方式

克隆仓库：

```bash
git clone https://github.com/Shawn-TV/Locus.git
cd Locus
```

创建虚拟环境：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

如果你的平台上 RDKit 用 pip 安装不顺，通常 Conda 更容易：

```bash
conda create -n locus python=3.11 rdkit pillow numpy opencv scikit-image -c conda-forge
conda activate locus
```

MolScribe 和 DECIMER 是可选项，只有需要用到时再按它们各自上游说明安装。模型权重请放在仓库外，不要提交进 Git。

## 准备数据

把比赛文件放到本地数据目录。常见结构如下：

```text
data/
├── pic/
│   ├── mol_0001.png
│   ├── mol_0002.png
│   └── ...
└── submission_template/
    └── submission.csv
```

`data/` 目录已经被 Git 忽略，因为比赛图片和模板受比赛主办方规则约束，不应该在这里重新分发。

创建工作目录：

```bash
mkdir -p work build
```

## 快速开始

先对原始图片调用 MolParser：

```bash
python tools/run_molparser_api.py \
  --input-dir data/pic \
  --output work/molparser_original.csv \
  --workers 6 \
  --resume
```

再对低置信度或疑难样本生成更安全的预处理版本并重新调用 MolParser：

```bash
python tools/run_molparser_variants.py \
  --input-dir data/pic \
  --from-csv work/molparser_original.csv \
  --threshold 0.95 \
  --include-aggressive \
  --output work/molparser_variants.csv
```

如果不用 MolScribe，可以先创建一个空的补充识别 CSV：

```bash
printf 'file_name,e_smiles,raw_smiles,confidence,error\n' > work/molscribe.csv
```

如果你已经安装 MolScribe 并准备好了 checkpoint，可以把它作为可选补充：

```bash
python tools/run_molscribe.py \
  --ckpt /path/to/molscribe.ckpt \
  --input-dir data/pic \
  --output work/molscribe.csv \
  --device cpu \
  --resume
```

生成一个初版提交包：

```bash
python tools/make_submission.py \
  --pic-dir data/pic \
  --template-csv data/submission_template/submission.csv \
  --molparser-csv work/molparser_original.csv work/molparser_variants.csv \
  --molscribe-csv work/molscribe.csv \
  --version V1 \
  --code-repo https://github.com/Shawn-TV/Locus \
  --out-dir build/V1 \
  --zip-path build/V1.zip
```

把预测的 E-SMILES 重新画出来，生成审阅表：

```bash
python tools/audit_v1_render_compare.py \
  --v1-zip build/V1.zip \
  --pic-dir data/pic \
  --out-dir reports/V1-audit \
  --molparser-csv work/molparser_original.csv work/molparser_variants.csv \
  --molscribe-csv work/molscribe.csv \
  --review-top 180
```

查看 `reports/V1-audit/review_sheets/` 里的图片。如果发现需要人工修正，可以参考 `examples/overrides.example.csv` 创建一个私人的 override CSV，然后重新生成：

```bash
python tools/make_submission.py \
  --pic-dir data/pic \
  --template-csv data/submission_template/submission.csv \
  --molparser-csv work/molparser_original.csv work/molparser_variants.csv \
  --molscribe-csv work/molscribe.csv \
  --overrides-csv /path/to/private_overrides.csv \
  --version V2 \
  --code-repo https://github.com/Shawn-TV/Locus \
  --out-dir build/V2 \
  --zip-path build/V2.zip
```

## 输出结果

最终 zip 里包含：

```text
submission.csv
meta.md
```

`submission.csv` 使用比赛要求的两列表头：

```csv
file_name,e_smiles
mol_0001.png,CCO<sep>
```

`meta.md` 会记录模型/工具使用情况、代码仓库链接和方法说明。

## 关于 E-SMILES

Locus 会尽量保留 MolParser 返回的 E-SMILES caption，因为比赛要求的是 Extended-SMILES，不只是普通 SMILES。合并脚本也会确保非空结果带有 `<sep>`，并保留原子标签等 E-SMILES 扩展标记。

缩写基团、Markush 标注、复杂桥环、立体化学和扫描图里的文字标签仍然可能需要人工复核。审计脚本就是为了这个循环设计的：回画、对比、检查、修正、重新打包。

## 第三方工具与许可证

Locus 使用 MIT License 开源。第三方工具仍然遵循它们自己的许可证和使用条款：

- MolParser public OCSR service：由 DP Technology 提供。
- MolScribe：可选的开源模型/checkpoint，遵循其上游项目许可证。
- RDKit：开源化学信息学工具包。
- DECIMER：可选 OCSR 辅助工具，遵循其上游项目许可证。
- Bohrium 比赛数据：由比赛主办方提供，本仓库不重新分发。

## 局限性

这个工具包可以帮助组织 OCSR 比赛流程，但不能保证所有化学结构都自动正确。低置信度识别、Markush 结构、立体化学、桥环、缩写基团和扫描标签仍然需要认真人工审阅。
