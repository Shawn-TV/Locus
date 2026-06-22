# Locus 复查工具使用说明

这份说明只讲怎么使用本地复查工具，不讲代码细节。

## 这个工具是做什么的

它用来人工检查和修改 `submission.csv` 里的 E-SMILES。

界面左边是图片列表，中间是原图，右边是当前 E-SMILES 重画出来的分子图。你改 E-SMILES 后，右边会自动重画，并且会自动保存到工作稿。

## 文件夹应该长什么样

正常情况下，项目里至少要有这些文件：

```text
Locus/
├── Start-Locus-Review.bat          Windows 双击启动
├── Start-Locus-Review.command      macOS 双击启动
├── start-locus-review.sh           Linux/终端启动
├── data/
│   └── pic/                        原始分子图片，里面是 mol_*.png
├── Track1-V1-final-submit/
│   ├── submission.csv              原始 V1 final 结果
│   └── meta.md
├── final_top300_review/
│   └── top300_risk_list.csv        前 300 张风险排序
└── final_review_editor/            自动生成，保存你的人工修改
```

如果启动时提示找不到图片，把 `mol_0001.png`、`mol_0002.png` 这类图片放到 `data/pic/` 里面。

## Windows 怎么启动

双击：

```text
Start-Locus-Review.bat
```

如果你拿到的是 `Locus-review-tool-Windows-portable.zip`，它已经自带 Python 运行环境，不需要自己安装 Python。解压后双击上面的 `.bat` 文件即可。不要关掉弹出的黑色窗口，浏览器会自动打开复查页面。

如果它提示找不到内置 runtime，说明你拿到的可能是源码版，不是 Windows 便携版。请换用：

```text
Locus-review-tool-Windows-portable.zip
```

## macOS 怎么启动

双击：

```text
Start-Locus-Review.command
```

如果系统提示不能打开，可以右键这个文件，然后选“打开”。

## Linux 怎么启动

在项目文件夹里运行：

```bash
./start-locus-review.sh
```

## 打开后怎么看

页面打开后，左侧会显示：

- `重点300`：风险最高的前 300 张，建议优先看。
- `全部`：全量 4000 张。
- `已修改`：你已经改过 E-SMILES 的图片。
- `存疑`：你手动标记为需要再看的图片。

注意：虽然人工重点看前 300 张，但最后导出的 `submission.csv` 仍然是全量 4000 行。

## 怎么修改

1. 在左边点一张图片。
2. 看中间原图。
3. 在右侧 `E-SMILES 编辑` 里修改文本。
4. 右边 RDKit 图会自动重画。
5. 停止输入后会自动保存，不需要每次手动点保存。

如果想立刻保存，也可以点 `保存工作稿`。

## 完成、存疑、备注是什么意思

- `完成这张`：表示你已经人工看过。
- `存疑`：表示这张以后还要再查。
- `备注`：写给自己的说明，比如“手性需再查”“右侧环系疑似不对”。

这些信息会自动保存到：

```text
final_review_editor/review_state.json
```

## 修改结果保存在哪里

自动保存的文件在：

```text
final_review_editor/working_submission.csv
```

这个文件永远是全量 4000 行。

只记录你改过哪些的文件在：

```text
final_review_editor/review_overrides.csv
```

## 怎么导出提交包

点页面右侧的：

```text
导出工作 zip
```

会生成并保存：

```text
final_review_editor/working_submission_package.zip
```

新版界面里也可以点：

```text
下载 CSV
下载 zip
```

这两个按钮会触发浏览器下载；同时文件仍然会保存在 `final_review_editor/` 目录里。

这个 zip 里面包含：

```text
submission.csv
meta.md
```

## 怎么确认有没有修改

看左上角：

```text
0 修改 · 0 存疑
```

如果显示 `5 修改`，说明当前工作稿里有 5 张和原始 V1 final 不一样。

## 怎么从 V1 final 重新开始

如果想清空所有人工修改，把 `Track1-V1-final-submit/submission.csv` 复制覆盖到：

```text
final_review_editor/working_submission.csv
```

然后重启工具即可。

注意：这样会丢掉当前人工修改，请确认后再做。
