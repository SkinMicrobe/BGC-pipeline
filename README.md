# BGC-Pipeline

![Bioinformatics](https://img.shields.io/badge/Bioinformatics-BGC%20Analysis-blue)
![Python](https://img.shields.io/badge/Python-3.8+-yellow)
![License](https://img.shields.io/badge/License-MIT-green)

BGC-Pipeline 是一个用于批量处理生物合成基因簇（Biosynthetic Gene Cluster, BGC）分析的完整自动化流程 Skill。它封装了从基因组注释、BGC 挖掘、去冗余到定量分析的全过程，旨在为研究人员提供一站式的分析解决方案。
<p align="center">
  <img src="/workflow.png" alt="BGC-Pipeline Workflow" width="750">
</p>


## ✨ Features (核心特性)

- 🧬 **双模式支持**：无缝兼容原核生物（细菌/古菌）和真核生物（如真菌）的分析流程。
- 🚀 **端到端自动化**：从组装好的 Fasta/MAG 文件输入，到最终的定量结果输出，一键运行。
- 🔄 **智能断点续传**：每个步骤均会自动生成 `.task.complete` 标记文件，意外中断后可无缝恢复，无需从头计算。
- 🧩 **高度模块化**：深度学习预测 (`deepBGC`) 和功能注释 (`eggNOG`) 均设为可选模块，方便灵活配置。
- 📊 **批量与定量分析**：支持多样本并行处理，并自动将测序数据（Fastq）映射到 BGC 上完成丰度定量。

---

## 🛠 Software & Versions (依赖软件与版本)

本流程依赖以下核心生物信息学工具。为保证流程的最佳兼容性与重现性，建议使用以下版本：

| 软件名称 | 版本要求 | 用途说明 | 所属 Conda 环境 |
| :--- | :--- | :--- | :--- |
| **antiSMASH** | `v8.0.4` | 核心 BGC 挖掘与预测 | `antismash` |
| **BiG-SCAPE** | `v2.0.3` | BGC 序列相似性网络分析与去冗余 | `bigmap` |
| **deepBGC** | `v0.1.31` | 基于深度学习的 BGC 预测 (可选) | `deepbgc` |
| **Bakta** | Latest | 原核生物基因组高质量注释 | `bakta` |
| **BiG-Map** | Latest | 丰度定量分析 (Family & Map) | `bigmap` |
| **metaEUK** | Latest | 真核生物高通量蛋白预测 | `meta` |
| **eggNOG-mapper** | Latest | BGC 蛋白深度功能注释 (可选) | `emapper` |

---

## 🗺️ Workflow (分析流)

### Prokaryote Mode (原核生物模式)
```text
MAGs (FASTA) → Bakta 注释 → antiSMASH 挖掘 
                             ├─→ [可选: deepBGC 预测] ─→ [可选: eggNOG 注释]
                             └─→ 提取 GBK → BiG-Map (去冗余) → BiG-Map (定量) → 合并汇总
```
## Eukaryote Mode (真核生物模式)
```text
MAGs (FASTA) → metaEUK 蛋白预测 → antiSMASH (真菌模式) → 提取 GBK → BiG-Map (去冗余)
```

---


## 📥 Downloads (下载与安装)
提供多种下载方式，请根据您的使用环境选择：

方法一：Git 克隆 (推荐)
适合需要随时同步更新到最新版本的用户：
```text
Bash
git clone [https://github.com/SkinMicrobe/BGC-pipeline.git](https://github.com/SkinMicrobe/BGC-pipeline.git)
cd BGC-pipeline
```
方法二：使用 Wget (适合 HPC/Linux 服务器)
如果您的服务器未配置 Git，可以直接通过 wget 下载压缩包：
```text
Bash
wget [https://github.com/SkinMicrobe/BGC-pipeline/archive/refs/heads/main.zip](https://github.com/SkinMicrobe/BGC-pipeline/archive/refs/heads/main.zip)
unzip main.zip
cd BGC-pipeline-main
```
方法三：直接下载 ZIP 压缩包
适合在 Windows/Mac 网页端操作的用户：

点击本页面右上角的 Code 绿色按钮，选择 Download ZIP。

解压下载的文件并进入该目录。

---


## ⚙️ Environment Setup (环境配置)
1. 准备 Conda 环境
需依次创建流程所需的独立环境（详见 envs/ 目录）：
```text
Bash
# 基础环境
conda env create -f envs/bakta.yaml
conda env create -f envs/antismash.yaml
conda env create -f envs/bigmap.yaml
conda env create -f envs/meta.yaml  # 仅真核模式需要

# 可选环境
conda env create -f envs/deepbgc.yaml
conda env create -f envs/emapper.yaml
```
2. 配置数据库路径
编辑 config/default_params.yaml 文件，将本地数据库路径替换为实际绝对路径：
```text
YAML
databases:
  bakta_db: "/path/to/bakta/db"
  pfam_hmm: "/path/to/Pfam-A.hmm/Pfam-A.hmm" # BiG-SCAPE 必需
  eggnog_db: "/path/to/eggnog_db_new"
  uniref50_db: "/path/to/uniref50_mmseqs"
```

---


## 🤖 交互使用指南 (Skill Integration & Usage)
本流程不仅是一个传统的命令行脚本，还被设计为可被大语言模型（如 Claude, Gemini）直接调用的 智能技能 (Skill)，并提供极度友好的交互体验。

💬 场景一：作为 AI Agent Skill 触发 (自然语言交互)
当您在集成了本 Skill 的 AI 助手环境中，直接使用自然语言提及相关关键词（如 BGC, antiSMASH, BiG-Map, 基因簇分析），即可唤醒该流程。

提示词 (Prompt) 示例：

🗣️ "我有一个真菌的 MAG 组装数据，在 /data/fungi_mag 目录下，帮我跑一下完整的 BGC 挖掘和去冗余，结果保存到 /data/results。"

🗣️ "使用原核模式分析 /data/bacteria/ 下的 contigs，我有对应的双端测序 fastq 文件在 /data/reads/，需要顺便运行 deepBGC 和 eggNOG 注释，帮我调用 BGC-pipeline 处理一下。"


⌨️ 场景二：命令行交互式运行 (CLI 问答)
如果您在终端中手动运行，无需记忆复杂的长串参数，系统会通过问答式引导您完成配置：
```text
Please select mode:
1) prokaryote (原核生物)
2) eukaryote (真核生物)
Choice: 1

Enter input fasta directory: /path/to/assembly/output
Enter output directory: /path/to/results
Number of threads: 32

Enter fastq directory (for quantification): /path/to/fastq
Run deepBGC? [y/n]: y
Run eggNOG? [y/n]: y
(注意：原核测序 Fastq 必须以 .fastq.gz 结尾；真核组装 MAG 支持 .fa.gz 或 .fa 格式。)
```

---


## 📂 Output Structure (输出目录结构)
所有分析结果将按照模块有序整理，主控日志与运行状态统一保存在 0_pipeline_info 中。

```text
OUT_DIR/
├── 0_pipeline_info/           # 运行配置、时间日志与断点标记
├── prok_results/              # [原核模式结果]
│   ├── 1_bakta_annotation/    
│   ├── 2_antismash_out/       
│   ├── 3_deepbgc_out/         # (可选)
│   ├── 4_deepbgc_extracted/   # (可选)
│   ├── 5_eggnog_mapper_out/   # (可选)
│   ├── 6_bigmap_family/       # 去冗余结果与提取的 GBK
│   ├── 7_bigmap_map_quant/    # 单样本丰度定量结果
│   └── 8_bigmap_summary/      # 跨样本定量的汇总表格
└── euk_results/               # [真核模式结果]
    ├── 1_metaeuk_annotation/
    ├── 2_antismash_out/
    ├── 3_gbk_extracted/
    └── 4_bigmap_family/
```

---


📄 License
This project is licensed under the MIT License.

✉️ Contact & Support
如有任何问题或优化建议，欢迎在 GitHub 提交 Issue。
