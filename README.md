# qPCR Auto Analysis & Plot

Python 脚本，用于自动化处理 qPCR 数据并绘图。

数据来源于 **Bio-Rad CFX Maestro** 软件，导出步骤：`Export → Custom Export`，格式为 CSV。

## 目录结构

```
├── README.md
├── config.yaml              # 默认配置文件
├── requirements.txt
├── data/                    # 原始输入 CSV
│   └── admin_*.csv
├── output/                  # 结果 CSV 输出
│   └── *_qPCR_result.csv
├── figures/                 # 生成的图表
│   └── *_qPCR_plot.png
├── src/                     # 核心模块
│   ├── config.py            # 配置加载
│   ├── io.py                # 文件读写
│   ├── preprocessing.py     # 数据清洗与分组
│   ├── analysis.py          # ΔCq → 2^(-ΔCq) → 归一化
│   ├── statistics.py        # t-test 与显著性
│   └── visualization.py     # 绘图
├── scripts/
│   └── run_analysis.py      # CLI 入口
└── notebooks/
    └── qPCR_analysis.ipynb  # 交互式 Notebook
```

## 快速开始

### 安装依赖

```bash
pip install -r requirements.txt
```

### 命令行使用

```bash
# 查看可用的 Target 基因
python scripts/run_analysis.py --input data/your_file.csv --list-targets

# 完整分析（使用默认配置）
python scripts/run_analysis.py \
    --input data/your_file.csv \
    --config config.yaml

# 指定参数
python scripts/run_analysis.py \
    --input data/your_file.csv \
    --ref-gene Rp49 \
    --control DMSO \
    --group-labels "PFOA001=PFOA 0.01mg/L,PFOA1=PFOA 1mg/L" \
    --output-dir output/ \
    --figures-dir figures/

# 仅导出 CSV，不绘图
python scripts/run_analysis.py --input data/your_file.csv --no-plot
```

### 命令行参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i, --input` | 输入的 Bio-Rad CSV 文件（必填） | — |
| `-c, --config` | YAML 配置文件路径 | `config.yaml` |
| `--ref-gene` | 内参基因名称 | Rp49 |
| `--control` | 对照组标签 | DMSO |
| `--separator` | 样本名中组名与重复号的分隔符 | `#` |
| `--group-labels` | 分组显示名映射 `raw=display,...` | — |
| `--skiprows` | CSV 元数据行数 | 19 |
| `-o, --output-dir` | CSV 输出目录 | `output/` |
| `-f, --figures-dir` | 图表输出目录 | `figures/` |
| `--no-plot` | 跳过绘图 | — |
| `--list-targets` | 列出文件中的 Target 基因后退出 | — |
| `-v, --verbose` | 详细日志 | — |

### Notebook 使用

在 Jupyter 中打开 `notebooks/qPCR_analysis.ipynb`，修改输入文件路径和分组标签后，依次运行所有 Cell。

## 分析流程

1. **读取数据** — 从 Bio-Rad CSV 提取 Target、Sample、Cq 列
2. **分组** — 按 Target 分组，按 Sample 排序
3. **计算 ΔCq** — Target_Cq − Reference_Cq
4. **计算 2^(-ΔCq)** — 表达量倍数
5. **归一化** — 除以对照组均值
6. **统计检验** — 各处理组 vs 对照组 t-test
7. **绘图** — 散点图 + 条形图 + 显著性标注

## 实验设计要求

- 样本命名格式：`组名#重复号`（如 `DMSO#1`、`PFOA001#2`）
- CSV 文件来自 Bio-Rad CFX Maestro 的 Custom Export
- 内参基因和目标基因都在 Target 列中

## 版本

- **v2.0** — 模块化重构：支持任意数量目的基因、自动分组、CLI + Notebook 双入口、完整错误处理
- **v1.x** — 原始版本（见仓库历史记录）
