# BGC Pipeline Skill

用于批量处理生物合成基因簇(BGC)分析的完整流程。

## 概述

本skill封装了从基因组注释到BGC定量分析的全流程，支持原核生物和真核生物两种模式。

### 支持的模式

#### 原核生物模式
```
组装FASTA → Bakta注释 → antiSMASH → [deepBGC] → [eggNOG] →
提取GBK → BiG-Map family去冗余 → BiG-Map map定量 → 合并结果
```

#### 真核生物模式
```
MAG FASTA → metaEUK注释 → antiSMASH（真菌）→ 提取GBK → BiG-Map family去冗余
```

## 目录结构

```
antismash-bgc/
├── SKILL.md                  # Skill配置文件
├── README.md                 # 本文件
├── bin/                      # 核心脚本
│   ├── run_pipeline.sh       # 主控脚本
│   ├── module_prok.sh        # 原核模块
│   ├── module_euk.sh         # 真核模块
│   └── utils/
│       ├── extract_gbk.py    # 提取GBK工具
│       ├── extract_deepbgc.py # 提取deepBGC蛋白
│       └── metaeuk_filter.py # metaEUK去重工具
├── config/                   # 配置文件
│   ├── default_params.yaml   # 默认参数（数据库路径等）
│   └── user_input.yaml       # 用户输入模板
└── envs/                     # Conda环境配置
    ├── antismash.yaml
    ├── bakta.yaml
    ├── bigmap.yaml
    ├── deepbgc.yaml          # 可选
    ├── emapper.yaml          # 可选
    └── meta.yaml
```

## 使用方法

### 1. 环境准备

首先需要创建所需的Conda环境：

```bash
# 创建基础环境
conda env create -f envs/bakta.yaml
conda env create -f envs/antismash.yaml
conda env create -f envs/bigmap.yaml
conda env create -f envs/meta.yaml  # 真核模式需要

# 创建可选环境
conda env create -f envs/deepbgc.yaml  # deepBGC（可选）
conda env create -f envs/emapper.yaml  # eggNOG（可选）

# 下载数据库（需要提前准备好）
# - Bakta数据库
# - antiSMASH数据库（conda安装时会下载）
# - Pfam数据库（BiG-SCAPE）
# - UniRef50数据库（metaEUK）
# - eggNOG数据库（可选）
```

### 2. 配置数据库路径

编辑 `config/default_params.yaml`，修改数据库路径为你的实际路径：

```yaml
databases:
  bakta_db: "/path/to/bakta/db"
  pfam_hmm: "/path/to/Pfam-A.hmm/Pfam-A.hmm"
  eggnog_db: "/path/to/eggnog_db_new"
  uniref50_db: "/path/to/uniref50_mmseqs"
```

### 3. 运行Pipeline

```bash
cd bin
bash run_pipeline.sh
```

按照提示输入：
- 模式选择（原核/真核）
- 输入目录
- 输出目录
- 线程数
- 原核模式需额外输入fastq目录
- 原核模式可选是否运行deepBGC和eggNOG

### 4. 断点续传

如果运行中断，重新运行会自动跳过已完成步骤。完成标记文件位于：
- 全局步骤：`0_pipeline_info/checkpoints/`
- 样本级别：各输出目录下的 `.task.complete`

## 输出目录结构

### 原核模式输出
```
OUT_DIR/
├── 0_pipeline_info/
│   ├── run_config.log
│   ├── run_time.log
│   ├── checkpoints/
│   └── logs/
└── prok_results/
    ├── 1_bakta_annotation/      # Bakta注释结果
    ├── 2_antismash_out/         # antiSMASH输出
    ├── 3_deepbgc_out/           # deepBGC结果（可选）
    ├── 4_deepbgc_extracted/     # 提取的蛋白序列（可选）
    ├── 5_eggnog_mapper_out/     # eggNOG注释（可选）
    ├── 6_bigmap_family/
    │   ├── aggregated_gbks/     # 提取的region gbk汇总
    │   └── bigscape_mash_results/  # 去冗余结果
    ├── 7_bigmap_map_quant/      # 定量结果
    └── 8_bigmap_summary/        # 合并的定量表格
```

### 真核模式输出
```
OUT_DIR/
├── 0_pipeline_info/
└── euk_results/
    ├── 1_metaeuk_annotation/    # metaEUK结果
    ├── 2_antismash_out/         # antiSMASH输出
    ├── 3_gbk_extracted/         # 提取的GBK汇总
    └── 4_bigmap_family/         # 去冗余结果
```

## 依赖软件

| 软件 | 用途 | Conda环境 | 必需 |
|------|------|-----------|------|
| Bakta | 原核基因组注释 | bakta | 是（原核） |
| antiSMASH | BGC挖掘 | antismash | 是 |
| BiG-Map | 去冗余和定量 | bigmap | 是 |
| metaEUK | 真核蛋白预测 | meta | 是（真核） |
| deepBGC | 深度学习BGC预测 | deepbgc | 否 |
| eggNOG-mapper | 功能注释 | emapper | 否 |

## 文件命名规范

### 原核模式输入
- **contigs**: `sample_dir/final.contigs.fa` 或 `*.fa`/`*.fasta`
- **fastq**: 任意命名，以 `.fastq.gz` 结尾
  - 自动识别R1/R2（支持：`R1`、`_1`、`read1`等模式）
  - 支持单端和双端数据
  - 单样本时直接复制结果，不进行合并操作

### 真核模式输入
- **MAG**: `input_dir/sample_name.fa.gz` 或 `.fa`

## 可选步骤说明

### deepBGC
使用深度学习预测BGC，可以补充antiSMASH的结果。
- 输入：原始contigs fasta文件
- 输出：每个样本一个GBK文件

### eggNOG mapper
对deepBGC提取的BGC蛋白进行功能注释。
- 输入：deepBGC提取的蛋白序列（.faa）
- 输出：eggNOG注释结果

## 常见问题

### Q: conda环境创建失败
A: 确保使用的是较新版本的conda，可以尝试使用mamba替代。

### Q: 数据库下载慢或失败
A: antiSMASH数据库会在首次运行时自动下载，建议使用国内镜像或提前下载。

### Q: BiG-Map运行报错
A: 确保Pfam数据库路径正确，且文件名为 `Pfam-A.hmm`。

### Q: fastq文件命名问题
A: 现在支持任意命名，只要以 `.fastq.gz` 结尾即可。会自动识别R1/R2。

### Q: 单样本定量结果
A: 如果只有一个样本，不会进行合并操作，直接复制结果文件。

## 版本历史

- v1.0: 初始版本，支持原核和真核两种模式
- v1.1: 添加deepBGC和eggNOG可选步骤，简化fastq文件处理

## 许可证

本skill基于MIT许可证开源。

## 联系方式

如有问题或建议，请在GitHub issue中提出。