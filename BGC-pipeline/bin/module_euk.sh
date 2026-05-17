#!/bin/bash

# ==========================================
# 真核生物BGC分析模块（真菌）
# ==========================================

set -euo pipefail

# 加载配置
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-}"
MODE="${MODE:-eukaryote}"
INPUT_DIR="${INPUT_DIR:-}"
THREADS="${THREADS:-16}"
TMP_DIR="${TMP_DIR:-/tmp/bgcpipeline_$$}"

# 从配置文件读取数据库路径
CONFIG_DIR="$(dirname "$SCRIPT_DIR")/config"
DEFAULT_PARAMS="$CONFIG_DIR/default_params.yaml"

# 解析数据库路径
UNIREF50_DB=$(grep -A1 "uniref50_db:" "$DEFAULT_PARAMS" 2>/dev/null | tail -1 | sed 's/.*: "\(.*\)".*/\1/' || echo "/mnt/cephfs/s2z1/db/uniref50/uniref50_mmseqs")

# 定义各步骤目录
METAEUK_OUT="$OUTPUT_DIR/euk_results/1_metaeuk_annotation"
ANTISMASH_OUT="$OUTPUT_DIR/euk_results/2_antismash_out"
EXTRACTED_GBKS="$OUTPUT_DIR/euk_results/3_gbk_extracted"
BIGMAP_FAMILY_OUT="$OUTPUT_DIR/euk_results/4_bigmap_family"

# 日志函数
log_info() { echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_warn() { echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2; }
log_error() { echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2; }

# 记录运行时间
record_time() {
    local step="$1"
    local action="$2"
    local time_file="$OUTPUT_DIR/0_pipeline_info/run_time.log"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $step - $action" >> "$time_file"
}

# 检查步骤是否完成
is_step_complete() {
    local marker="$1"
    [[ "${SKIP_EXISTING:-true}" == "true" && -f "$marker" ]]
}

# 设置断点标记
mark_complete() {
    local marker_file="$1"
    touch "$marker_file"
}

# ==========================================
# 步骤1: metaEUK注释
# ==========================================
step1_metaeuk() {
    local step_name="step1_metaeuk"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤1: metaEUK注释已完成，跳过"
        return 0
    fi

    record_time "步骤1_metaEUK注释" "开始"
    log_info "开始步骤1: metaEUK蛋白预测"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate meta

    mkdir -p "$TMP_DIR"
    export TMPDIR="$TMP_DIR"

    shopt -s nullglob
    local processed=0
    local skipped=0

    for fasta in "$INPUT_DIR"/*.fa.gz "$INPUT_DIR"/*.fa; do
        if [[ ! -f "$fasta" ]]; then
            continue
        fi

        fname=$(basename "$fasta")
        # 去除扩展名
        name="${fname%.fa.gz}"
        name="${name%.fa}"

        outdir="$METAEUK_OUT/$name"
        marker_file="$outdir/.task.complete"

        if is_step_complete "$marker_file"; then
            log_info "样本 ${name} 已完成，跳过"
            ((skipped++))
            continue
        fi

        log_info "正在处理样本: ${name}"
        mkdir -p "$outdir"

        # 独立tmp目录
        sample_tmp="$TMP_DIR/tmp_${name}"
        mkdir -p "$sample_tmp"

        out_prefix="$outdir/predicted_proteins"

        if metaeuk easy-predict \
            "$fasta" \
            "$UNIREF50_DB" \
            "$out_prefix" \
            "$sample_tmp" \
            --threads "$THREADS" >> "$log_file" 2>&1; then
            mark_complete "$marker_file"
            ((processed++))
            log_info "样本 ${name} 处理完成"
            rm -rf "$sample_tmp"
        else
            log_error "样本 ${name} 处理失败"
            return 1
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤1_metaEUK注释" "完成"
    log_info "步骤1完成: 处理 ${processed} 个样本，跳过 ${skipped} 个样本"
}

# ==========================================
# 步骤2: antiSMASH（真菌模式）
# ==========================================
step2_antismash_euk() {
    local step_name="step2_antismash_euk"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤2: antiSMASH（真菌）已完成，跳过"
        return 0
    fi

    record_time "步骤2_antiSMASH_真菌" "开始"
    log_info "开始步骤2: antiSMASH挖掘（真菌模式）"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate antismash

    shopt -s nullglob
    local processed=0
    local skipped=0

    # 遍历metaEUK输出目录
    for MAG_DIR in "$METAEUK_OUT"/*/; do
        if [[ ! -d "$MAG_DIR" ]]; then
            continue
        fi

        MAG_NAME=$(basename "$MAG_DIR")

        # 原始GFF
        GFF="$MAG_DIR/predicted_proteins.gff"

        # 去重后的CDS GFF
        CDS_GFF="$MAG_DIR/predicted_proteins_CDS.gff"

        # 对应的fasta文件（从输入目录查找）
        FA=$(find "$INPUT_DIR" -name "${MAG_NAME}.fa.gz" -o -name "${MAG_NAME}.fa" 2>/dev/null | head -1)

        if [[ ! -f "$GFF" ]]; then
            log_warn "样本 ${MAG_NAME} 未找到GFF，跳过"
            continue
        fi

        if [[ ! -f "$FA" ]]; then
            log_warn "样本 ${MAG_NAME} 未找到fasta文件，跳过"
            continue
        fi

        # 检查是否已完成
        local sample_checkpoint="$ANTISMASH_OUT/${MAG_NAME}/.task.complete"
        if is_step_complete "$sample_checkpoint"; then
            log_info "样本 ${MAG_NAME} 已完成，跳过"
            ((skipped++))
            continue
        fi

        log_info "正在处理样本: ${MAG_NAME}"
        mkdir -p "$ANTISMASH_OUT/${MAG_NAME}"

        # 1. 生成去重CDS文件
        log_info "生成去重CDS GFF..."
        awk 'BEGIN{FS=OFS="\t"}
             /^#/ {print; next}
             $3=="CDS" {
                 key = $1 FS $4 FS $5 FS $7
                 if (!seen[key]++) print
             }' "$GFF" > "$CDS_GFF"

        # 2. 运行antiSMASH
        log_info "运行antiSMASH..."
        if antismash "$FA" \
            --taxon fungi \
            --genefinding-tool none \
            --genefinding-gff3 "$CDS_GFF" \
            --output-dir "$ANTISMASH_OUT/${MAG_NAME}" \
            --cpus "$THREADS" >> "$log_file" 2>&1; then
            mark_complete "$sample_checkpoint"
            ((processed++))
            log_info "样本 ${MAG_NAME} 处理完成"
        else
            log_error "样本 ${MAG_NAME} 处理失败"
            return 1
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤2_antiSMASH_真菌" "完成"
    log_info "步骤2完成: 处理 ${processed} 个样本，跳过 ${skipped} 个样本"
}

# ==========================================
# 步骤3: 提取GBK文件
# ==========================================
step3_extract_gbk() {
    local step_name="step3_extract_gbk"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"

    if is_step_complete "$checkpoint"; then
        log_info "步骤3: 提取GBK已完成，跳过"
        return 0
    fi

    record_time "步骤3_提取GBK" "开始"
    log_info "开始步骤3: 提取region gbk文件"

    mkdir -p "$EXTRACTED_GBKS"

    shopt -s nullglob
    local total_files=0

    for SAMPLE_DIR in "$ANTISMASH_OUT"/*/; do
        SAMPLE_NAME=$(basename "${SAMPLE_DIR%/}")

        for GBK_FILE in "$SAMPLE_DIR"/*.region*.gbk; do
            if [[ -f "$GBK_FILE" ]]; then
                ORIGINAL_NAME=$(basename "$GBK_FILE")
                NEW_NAME="${SAMPLE_NAME}_${ORIGINAL_NAME}"
                cp "$GBK_FILE" "${EXTRACTED_GBKS}/${NEW_NAME}"
                ((total_files++))
            fi
        done
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤3_提取GBK" "完成"
    log_info "步骤3完成: 共提取 ${total_files} 个region gbk文件"
}

# ==========================================
# 步骤4: BiG-Map family（去冗余）
# ==========================================
step4_bigmap_family() {
    local step_name="step4_bigmap_family"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤4: BiG-Map family已完成，跳过"
        return 0
    fi

    record_time "步骤4_BiG-Map_family" "开始"
    log_info "开始步骤4: BiG-Map family去冗余"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate bigmap

    # 获取Pfam路径
    local pfam_hmm=$(grep -A1 "pfam_hmm:" "$DEFAULT_PARAMS" 2>/dev/null | tail -1 | sed 's/.*: "\(.*\)".*/\1/' || echo "")

    if [[ ! -f "$pfam_hmm" ]]; then
        log_error "未找到Pfam数据库: $pfam_hmm"
        return 1
    fi

    mkdir -p "$BIGMAP_FAMILY_OUT"

    # 获取bigmap family线程数（默认48）
    local bigmap_threads=$(grep -A5 "bigmap_family:" "$DEFAULT_PARAMS" | grep -oP '(?<=: )\d+' || echo "48")

    if bigmap-family \
        -D "$EXTRACTED_GBKS" \
        -b bigscape \
        -pf "$pfam_hmm" \
        -p "$bigmap_threads" \
        -O "$BIGMAP_FAMILY_OUT" >> "$log_file" 2>&1; then
        mark_complete "$checkpoint"
        record_time "步骤4_BiG-Map_family" "完成"
        log_info "步骤4完成"
    else
        log_error "BiG-Map family运行失败"
        return 1
    fi
}

# ==========================================
# 主函数
# ==========================================
main() {
    log_info "=========================================="
    log_info "   真核生物BGC分析模块启动"
    log_info "=========================================="

    step1_metaeuk || exit 1
    step2_antismash_euk || exit 1
    step3_extract_gbk || exit 1
    step4_bigmap_family || exit 1

    log_info "=========================================="
    log_info "   真核生物BGC分析模块完成"
    log_info "=========================================="
}

main "$@"