---
name: BGC-pipeline
description: 批量处理生物合成基因簇(BGC)分析全流程。当用户需要进行antiSMASH BGC挖掘、BGC去冗余、定量分析时使用此skill。支持原核生物(Bakta→antiSMASH→[deepBGC]→[eggNOG]→BiG-Map family→BiG-Map map)和真核生物(metaEUK→antiSMASH→BiG-Map family)两种模式。触发条件：用户提到BGC、biosynthetic gene cluster、antiSMASH、BiG-Map、基因簇分析、次级代谢产物等。
---

# BGC Pipeline Skill

用于批量处理生物合成基因簇(BGC)分析的完整流程。

## 支持的模式

### 原核生物模式
完整流程包括：
1. **Bakta注释** - 对组装contigs进行基因功能注释
2. **antiSMASH挖掘** - 使用Bakta注释结果进行BGC预测
3. **deepBGC预测（可选）** - 使用深度学习预测BGC
4. **提取deepBGC序列（可选）** - 从deepBGC输出提取蛋白序列
5. **eggNOG注释（可选）** - 对BGC蛋白进行功能注释
6. **提取GBK文件** - 汇总antiSMASH输出的region gbk文件
7. **BiG-Map family去冗余** - 使用Mash和BiG-SCAPE 2.0进行去冗余
8. **BiG-Map map定量** - 将原始测序数据映射到BGC上进行定量
9. **合并结果** - 汇总所有样本的定量表格（单样本直接复制）

### 真核生物模式
完整流程包括：
1. **metaEUK注释** - 对MAG进行蛋白预测
2. **antiSMASH挖掘** - 去重CDS后进行真菌BGC预测
3. **提取GBK文件** - 汇总所有region gbk文件
4. **BiG-Map family去冗余** - 去冗余分析（无定量步骤）

## 执行流程

### 第一步：收集用户输入

运行主控脚本 `bin/run_pipeline.sh`，它会提示用户输入：

**必填参数：**
- 模式选择：`prokaryote`（原核）或 `eukaryote`（真核）
- 输入fasta/MAG目录路径
- 输出目录路径
- 线程数

**原核模式额外参数：**
- 测序fastq文件目录路径（用于定量，支持单端和双端，以.fastq.gz结尾）
- 是否运行deepBGC（可选）
- 是否运行eggNOG（可选）

**真核模式无需fastq输入**

### 第二步：配置检查与确认

脚本会：
1. 检查各软件Conda环境是否可激活
2. 检查数据库路径是否存在
3. 创建输出目录结构
4. 保存运行配置到 `0_pipeline_info/run_config.log`

### 第三步：按顺序执行各步骤

每个步骤完成后会生成 `.task.complete` 标记文件，支持断点续传。

## 输出目录结构

```
OUT_DIR/
├── 0_pipeline_info/           # 全局配置与运行信息
│   ├── run_config.log         # 运行配置快照
│   ├── run_time.log           # 每步起止时间
│   ├── checkpoints/           # 断点续传标记
│   └── logs/                  # 详细日志
│
├── prok_results/              # 原核生物结果目录
│   ├── 1_bakta_annotation/    # Bakta注释结果（含.gbff）
│   ├── 2_antismash_out/       # antiSMASH输出
│   ├── 3_deepbgc_out/         # [可选] deepBGC预测结果
│   ├── 4_deepbgc_extracted/   # [可选] 提取的BGC蛋白序列
│   ├── 5_eggnog_mapper_out/   # [可选] eggNOG注释结果
│   ├── 6_bigmap_family/
│   │   ├── aggregated_gbks/   # 提取的region gbk汇总
│   │   └── bigscape_mash_results/  # 去冗余结果
│   ├── 7_bigmap_map_quant/    # 定量结果
│   └── 8_bigmap_summary/      # 合并的定量表格
│
└── euk_results/               # 真核生物结果目录
    ├── 1_metaeuk_annotation/  # metaEUK蛋白预测结果
    ├── 2_antismash_out/       # antiSMASH输出
    ├── 3_gbk_extracted/       # 提取的GBK汇总
    └── 4_bigmap_family/       # 去冗余结果
```

## 配置文件

### config/default_params.yaml
默认参数配置，包括：
- 各软件Conda环境名
- 数据库绝对路径（Bakta、Pfam、UniRef50、eggNOG）
- 默认线程数

### config/user_input.yaml
用户每次运行前填写，包含：
- 输入输出路径
- 模式选择
- CPU设置
- 可选步骤开关

## Conda环境

需要在以下环境中安装对应软件：
- `bakta` - Bakta注释工具
- `antismash` - antiSMASH BGC挖掘
- `bigmap` - BiG-Map family和BiG-Map map
- `meta` - metaEUK（仅真核模式）
- `deepbgc` - deepBGC（原核可选）
- `emapper` - eggNOG mapper（原核可选）

## 工具脚本

### bin/utils/extract_gbk.py
提取antiSMASH输出的 `.region*.gbk` 文件并重命名（添加样本前缀防覆盖）

### bin/utils/extract_deepbgc.py
从deepBGC输出的GBK文件中提取BGC区域的CDS蛋白序列

### bin/utils/metaeuk_filter.py
处理metaEUK结果，去重CDS并生成GFF3

## 断点续传

每个步骤完成后生成对应的 `.task.complete` 标记文件。重新运行时会检查标记文件，跳过已完成步骤。

## 使用示例

```bash
# 激活主控脚本
cd bin
bash run_pipeline.sh

# 按提示输入参数
# 模式: prokaryote
# 输入fasta目录: /path/to/megahit/output
# 输出目录: /path/to/results
# 线程数: 32
# fastq目录: /path/to/fastq
# 是否运行deepBGC: n
# 是否运行eggNOG: n
```

## 注意事项

1. 确保所有conda环境都已正确配置
2. 确保数据库路径在 `config/default_params.yaml` 中正确设置
3. 原核模式fastq文件需以 `.fastq.gz` 结尾，支持单端和双端数据
4. 真核模式需要MAG文件格式为 `.fa.gz` 或 `.fa`
5. deepBGC和eggNOG为可选步骤，不影响主流程
6. 大样本建议使用高线程数配置