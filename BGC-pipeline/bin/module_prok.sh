#!/bin/bash

# ==========================================
# 原核生物BGC分析模块
# ==========================================

set -euo pipefail

# 加载配置
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-}"
MODE="${MODE:-prokaryote}"
INPUT_DIR="${INPUT_DIR:-}"
FASTQ_DIR="${FASTQ_DIR:-}"
THREADS="${THREADS:-32}"
TMP_DIR="${TMP_DIR:-/tmp/bgcpipeline_$$}"
RUN_DEEPBGC="${RUN_DEEPBGC:-false}"
RUN_EGGNOG="${RUN_EGGNOG:-false}"

# 从配置文件读取数据库路径
CONFIG_DIR="$(dirname "$SCRIPT_DIR")/config"
DEFAULT_PARAMS="$CONFIG_DIR/default_params.yaml"

# 解析数据库路径
BAKTA_DB=$(grep "bakta_db:" "$DEFAULT_PARAMS" | sed 's/.*: "\(.*\)".*/\1/' || echo "/mnt/cephfs/s2z1/db/bakta/full/db")
PFAM_HMM=$(grep "pfam_hmm:" "$DEFAULT_PARAMS" | sed 's/.*: "\(.*\)".*/\1/' || echo "")
EGGNOG_DB=$(grep "eggnog_db:" "$DEFAULT_PARAMS" | sed 's/.*: "\(.*\)".*/\1/' || echo "/mnt/cephfs/s2z1/db/eggnog_db_new")

# 定义各步骤目录
BAKTA_OUT="$OUTPUT_DIR/prok_results/1_bakta_annotation"
ANTISMASH_OUT="$OUTPUT_DIR/prok_results/2_antismash_out"
DEEPBGC_OUT="$OUTPUT_DIR/prok_results/3_deepbgc_out"
DEEPBGC_EXTRACTED="$OUTPUT_DIR/prok_results/4_deepbgc_extracted"
EGGNOG_OUT="$OUTPUT_DIR/prok_results/5_eggnog_mapper_out"
AGGREGATED_GBKS="$OUTPUT_DIR/prok_results/6_bigmap_family/aggregated_gbks"
BIGMAP_FAMILY_OUT="$OUTPUT_DIR/prok_results/6_bigmap_family/bigscape_mash_results"
BIGMAP_MAP_OUT="$OUTPUT_DIR/prok_results/7_bigmap_map_quant"
BIGMAP_SUMMARY="$OUTPUT_DIR/prok_results/8_bigmap_summary"

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
# 步骤1: Bakta注释
# ==========================================
step1_bakta() {
    local step_name="step1_bakta"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤1: Bakta注释已完成，跳过"
        return 0
    fi

    record_time "步骤1_Bakta注释" "开始"
    log_info "开始步骤1: Bakta注释"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate bakta

    # 开启nullglob
    shopt -s nullglob

    local processed=0
    local skipped=0

    # 只处理final.contigs.fa（跳过其他.fa文件）
    local fasta_files=()
    if [[ -f "$INPUT_DIR/final.contigs.fa" ]]; then
        fasta_files+=("$INPUT_DIR/final.contigs.fa")
    else
        for dir in "$INPUT_DIR"/*/; do
            if [[ -f "$dir/final.contigs.fa" ]]; then
                fasta_files+=("$dir/final.contigs.fa")
            fi
        done
    fi

    for FASTA in "${fasta_files[@]}"; do
        # 提取样本名
        if [[ "$FASTA" == */* ]]; then
            SAMPLE_NAME=$(basename "$(dirname "$FASTA")")
        else
            SAMPLE_NAME=$(basename "$FASTA" .fa)
            SAMPLE_NAME=$(basename "$SAMPLE_NAME" .fasta)
        fi

        SAMPLE_OUT_DIR="${BAKTA_OUT}/${SAMPLE_NAME}"
        local sample_checkpoint="${SAMPLE_OUT_DIR}/.task.complete"

        if is_step_complete "$sample_checkpoint"; then
            log_info "样本 ${SAMPLE_NAME} 已完成，跳过"
            ((skipped++))
            continue
        fi

        log_info "正在处理样本: ${SAMPLE_NAME}"

        if bakta \
            --db "$BAKTA_DB" \
            --output "$SAMPLE_OUT_DIR" \
            --prefix "$SAMPLE_NAME" \
            --threads "$THREADS" \
            --force \
            "$FASTA" >> "$log_file" 2>&1; then
            mark_complete "$sample_checkpoint"
            ((processed++))
            log_info "样本 ${SAMPLE_NAME} 处理完成"
        else
            log_error "样本 ${SAMPLE_NAME} 处理失败"
            return 1
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤1_Bakta注释" "完成"
    log_info "步骤1完成: 处理 ${processed} 个样本，跳过 ${skipped} 个样本"
}

# ==========================================
# 步骤2: antiSMASH（细菌模式）
# ==========================================
step2_antismash() {
    local step_name="step2_antismash"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤2: antiSMASH已完成，跳过"
        return 0
    fi

    record_time "步骤2_antiSMASH" "开始"
    log_info "开始步骤2: antiSMASH挖掘（细菌模式）"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate antismash

    shopt -s nullglob

    local processed=0
    local skipped=0

    for GBFF in "$BAKTA_OUT"/*/*.gbff; do
        if [[ ! -f "$GBFF" ]]; then
            continue
        fi

        SAMPLE_NAME=$(basename "$(dirname "$GBFF")")
        SAMPLE_OUT_DIR="${ANTISMASH_OUT}/${SAMPLE_NAME}"
        local sample_checkpoint="${SAMPLE_OUT_DIR}/.task.complete"

        if is_step_complete "$sample_checkpoint"; then
            log_info "样本 ${SAMPLE_NAME} 已完成，跳过"
            ((skipped++))
            continue
        fi

        log_info "正在处理样本: ${SAMPLE_NAME}"
        mkdir -p "$SAMPLE_OUT_DIR"

        if antismash \
            --cpus "$THREADS" \
            --taxon bacteria \
            --output-dir "$SAMPLE_OUT_DIR" \
            --genefinding-tool none \
            --cb-general --cb-knownclusters --cb-subclusters \
            --asf --pfam2go \
            "$GBFF" >> "$log_file" 2>&1; then
            mark_complete "$sample_checkpoint"
            ((processed++))
            log_info "样本 ${SAMPLE_NAME} 处理完成"
        else
            log_error "样本 ${SAMPLE_NAME} 处理失败"
            return 1
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤2_antiSMASH" "完成"
    log_info "步骤2完成: 处理 ${processed} 个样本，跳过 ${skipped} 个样本"
}

# ==========================================
# 步骤3: deepBGC（可选）
# ==========================================
step3_deepbgc() {
    if [[ "$RUN_DEEPBGC" != "true" ]]; then
        log_info "步骤3: deepBGC未启用，跳过"
        return 0
    fi

    local step_name="step3_deepbgc"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤3: deepBGC已完成，跳过"
        return 0
    fi

    record_time "步骤3_deepBGC" "开始"
    log_info "开始步骤3: deepBGC BGC预测"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate deepbgc

    shopt -s nullglob
    local processed=0
    local skipped=0

    # 输入原始fasta文件
    for FASTA in "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fasta; do
        if [[ ! -f "$FASTA" ]]; then
            continue
        fi

        filename=$(basename "$FASTA" .fa)
        filename=$(basename "$filename" .fasta)

        output_dir="$DEEPBGC_OUT/$filename"
        local sample_checkpoint="${output_dir}/task.complete"

        if is_step_complete "$sample_checkpoint"; then
            log_info "样本 ${filename} 已完成，跳过"
            ((skipped++))
            continue
        fi

        log_info "正在处理样本: ${filename}"
        mkdir -p "$output_dir"

        if deepbgc pipeline \
            --prodigal-meta-mode \
            --output "$output_dir" \
            --detector deepbgc \
            --classifier product_class \
            --classifier product_activity \
            --score 0.5 \
            --minimal-output \
            "$FASTA" >> "$log_file" 2>&1; then
            mark_complete "$sample_checkpoint"
            ((processed++))
            log_info "样本 ${filename} 处理完成"
        else
            log_warn "样本 ${filename} 处理失败"
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤3_deepBGC" "完成"
    log_info "步骤3完成: 处理 ${processed} 个样本，跳过 ${skipped} 个样本"
}

# ==========================================
# 步骤4: 提取deepBGC序列（可选）
# ==========================================
step4_extract_deepbgc() {
    if [[ "$RUN_DEEPBGC" != "true" ]]; then
        log_info "步骤4: 提取deepBGC序列未启用，跳过"
        return 0
    fi

    local step_name="step4_extract_deepbgc"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤4: 提取deepBGC序列已完成，跳过"
        return 0
    fi

    record_time "步骤4_提取deepBGC序列" "开始"
    log_info "开始步骤4: 提取deepBGC蛋白序列"

    mkdir -p "$DEEPBGC_EXTRACTED"

    shopt -s nullglob
    local processed=0
    local total_proteins=0

    for sample_dir in "$DEEPBGC_OUT"/*/; do
        if [[ ! -d "$sample_dir" ]]; then
            continue
        fi

        sample_name=$(basename "$sample_dir")
        local sample_checkpoint="${DEEPBGC_EXTRACTED}/${sample_name}.complete"

        if is_step_complete "$sample_checkpoint"; then
            log_info "样本 ${sample_name} 已完成，跳过"
            ((processed++))
            continue
        fi

        # 查找GBK文件
        gbk_file=$(find "$sample_dir" -name "*.gbk" -type f | head -1)

        if [[ ! -f "$gbk_file" ]]; then
            log_warn "样本 ${sample_name} 未找到GBK文件，跳过"
            continue
        fi

        log_info "正在处理样本: ${sample_name}"
        output_faa="${DEEPBGC_EXTRACTED}/${sample_name}_bgc_proteins.faa"

        # 运行提取脚本
        if python3 "$SCRIPT_DIR/utils/extract_deepbgc.py" \
            "$gbk_file" "$output_faa" >> "$log_file" 2>&1; then
            mark_complete "$sample_checkpoint"
            ((processed++))
            log_info "样本 ${sample_name} 处理完成"
        else
            log_warn "样本 ${sample_name} 提取失败"
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤4_提取deepBGC序列" "完成"
    log_info "步骤4完成: 处理 ${processed} 个样本"
}

# ==========================================
# 步骤5: eggNOG mapper（可选）
# ==========================================
step5_eggnog_mapper() {
    if [[ "$RUN_EGGNOG" != "true" ]]; then
        log_info "步骤5: eggNOG mapper未启用，跳过"
        return 0
    fi

    local step_name="step5_eggnog_mapper"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤5: eggNOG mapper已完成，跳过"
        return 0
    fi

    record_time "步骤5_eggNOG_mapper" "开始"
    log_info "开始步骤5: eggNOG mapper功能注释"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate emapper

    mkdir -p "$EGGNOG_OUT"

    shopt -s nullglob
    local processed=0
    local skipped=0

    for faa_file in "$DEEPBGC_EXTRACTED"/*_bgc_proteins.faa; do
        if [[ ! -f "$faa_file" ]]; then
            continue
        fi

        faa_basename=$(basename "$faa_file" _bgc_proteins.faa)
        sample_output_dir="${EGGNOG_OUT}/${faa_basename}"
        local sample_checkpoint="${sample_output_dir}/.task.complete"

        if is_step_complete "$sample_checkpoint"; then
            log_info "样本 ${faa_basename} 已完成，跳过"
            ((skipped++))
            continue
        fi

        log_info "正在处理样本: ${faa_basename}"
        mkdir -p "$sample_output_dir"

        if emapper.py \
            -m diamond \
            -i "$faa_file" \
            -o "${faa_basename}_bgc" \
            --data_dir "$EGGNOG_DB" \
            --output_dir "$sample_output_dir" \
            --cpu "$THREADS" \
            --itype proteins \
            --override >> "$log_file" 2>&1; then
            mark_complete "$sample_checkpoint"
            ((processed++))
            log_info "样本 ${faa_basename} 处理完成"
        else
            log_warn "样本 ${faa_basename} 处理失败"
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤5_eggNOG_mapper" "完成"
    log_info "步骤5完成: 处理 ${processed} 个样本，跳过 ${skipped} 个样本"
}

# ==========================================
# 步骤6: 提取GBK文件
# ==========================================
step6_extract_gbk() {
    local step_name="step6_extract_gbk"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"

    if is_step_complete "$checkpoint"; then
        log_info "步骤6: 提取GBK已完成，跳过"
        return 0
    fi

    record_time "步骤6_提取GBK" "开始"
    log_info "开始步骤6: 提取region gbk文件"

    mkdir -p "$AGGREGATED_GBKS"

    shopt -s nullglob
    local total_files=0

    for SAMPLE_DIR in "$ANTISMASH_OUT"/*/; do
        SAMPLE_NAME=$(basename "${SAMPLE_DIR%/}")

        for GBK_FILE in "$SAMPLE_DIR"/*.region*.gbk; do
            if [[ -f "$GBK_FILE" ]]; then
                ORIGINAL_NAME=$(basename "$GBK_FILE")
                NEW_NAME="${SAMPLE_NAME}_${ORIGINAL_NAME}"
                cp "$GBK_FILE" "${AGGREGATED_GBKS}/${NEW_NAME}"
                ((total_files++))
            fi
        done
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤6_提取GBK" "完成"
    log_info "步骤6完成: 共提取 ${total_files} 个region gbk文件"
}

# ==========================================
# 步骤7: BiG-Map family（去冗余）
# ==========================================
step7_bigmap_family() {
    local step_name="step7_bigmap_family"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤7: BiG-Map family已完成，跳过"
        return 0
    fi

    record_time "步骤7_BiG-Map_family" "开始"
    log_info "开始步骤7: BiG-Map family去冗余"

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate bigmap

    # 获取bigmap family线程数（默认48）
    local bigmap_threads=$(grep -A5 "bigmap_family:" "$DEFAULT_PARAMS" | grep -oP '(?<=: )\d+' || echo "48")

    mkdir -p "$BIGMAP_FAMILY_OUT"

    if bigmap-family \
        -D "$AGGREGATED_GBKS" \
        -b bigscape \
        -pf "$PFAM_HMM" \
        -p "$bigmap_threads" \
        -O "$BIGMAP_FAMILY_OUT" >> "$log_file" 2>&1; then
        mark_complete "$checkpoint"
        record_time "步骤7_BiG-Map_family" "完成"
        log_info "步骤7完成"
    else
        log_error "BiG-Map family运行失败"
        return 1
    fi
}

# ==========================================
# 步骤8: BiG-Map map（定量）
# ==========================================
step8_bigmap_map() {
    local step_name="step8_bigmap_map"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"
    local log_file="$OUTPUT_DIR/0_pipeline_info/logs/${step_name}.log"

    if is_step_complete "$checkpoint"; then
        log_info "步骤8: BiG-Map map已完成，跳过"
        return 0
    fi

    record_time "步骤8_BiG-Map_map" "开始"
    log_info "开始步骤8: BiG-Map map定量"

    # 检查参考文件
    if [[ ! -f "$BIGMAP_FAMILY_OUT/BiG-MAP.GCF.json" ]]; then
        log_error "未找到 BiG-MAP.GCF.json，请确认步骤7已完成"
        return 1
    fi

    # 激活环境
    eval "$(conda shell.bash hook)"
    conda activate bigmap

    # 获取bigmap map线程数（默认36）
    local bigmap_map_threads=$(grep -A5 "bigmap_map:" "$DEFAULT_PARAMS" | grep -oP '(?<=: )\d+' || echo "36")

    mkdir -p "$BIGMAP_MAP_OUT"

    shopt -s nullglob
    local processed=0
    local skipped=0

    # 遍历fastq目录
    for fastq_dir in "$FASTQ_DIR"/*/; do
        if [[ ! -d "$fastq_dir" ]]; then
            continue
        fi

        sample_name=$(basename "$fastq_dir")
        out_dir="$BIGMAP_MAP_OUT/${sample_name}"
        local sample_checkpoint="${out_dir}/.task.complete"

        if is_step_complete "$sample_checkpoint"; then
            log_info "样本 ${sample_name} 已完成，跳过"
            ((skipped++))
            continue
        fi

        # 查找fastq文件（支持单端和双端）
        r1_files=($(find "$fastq_dir" -maxdepth 1 -name "*.fastq.gz" | grep -E "R1|_1|read1" | head -1))
        r2_files=($(find "$fastq_dir" -maxdepth 1 -name "*.fastq.gz" | grep -E "R2|_2|read2" | head -1))

        mate1_file=""
        mate2_file=""

        if [[ ${#r1_files[@]} -gt 0 ]]; then
            mate1_file="${r1_files[0]}"
        fi

        if [[ ${#r2_files[@]} -gt 0 ]]; then
            mate2_file="${r2_files[0]}"
        fi

        # 如果没有找到R1/R2，尝试找任意两个fastq.gz文件
        if [[ -z "$mate1_file" ]]; then
            all_files=($(find "$fastq_dir" -maxdepth 1 -name "*.fastq.gz" | sort))
            if [[ ${#all_files[@]} -ge 2 ]]; then
                mate1_file="${all_files[0]}"
                mate2_file="${all_files[1]}"
            elif [[ ${#all_files[@]} -eq 1 ]]; then
                mate1_file="${all_files[0]}"
                log_warn "样本 ${sample_name} 只有单端数据"
            fi
        fi

        if [[ -z "$mate1_file" ]]; then
            log_warn "样本 ${sample_name} 未找到fastq文件，跳过"
            continue
        fi

        log_info "正在处理样本: ${sample_name}"
        mkdir -p "$out_dir"

        if [[ -n "$mate2_file" ]]; then
            # 双端数据
            if bigmap-map \
                -I1 "$mate1_file" \
                -I2 "$mate2_file" \
                -F "$BIGMAP_FAMILY_OUT" \
                -th "$bigmap_map_threads" \
                -s sensitive-local \
                -O "$out_dir" >> "$log_file" 2>&1; then
                mark_complete "$sample_checkpoint"
                ((processed++))
                log_info "样本 ${sample_name} 处理完成（双端）"
            else
                log_warn "样本 ${sample_name} 处理失败"
            fi
        else
            # 单端数据
            if bigmap-map \
                -I1 "$mate1_file" \
                -F "$BIGMAP_FAMILY_OUT" \
                -th "$bigmap_map_threads" \
                -s sensitive-local \
                -O "$out_dir" >> "$log_file" 2>&1; then
                mark_complete "$sample_checkpoint"
                ((processed++))
                log_info "样本 ${sample_name} 处理完成（单端）"
            else
                log_warn "样本 ${sample_name} 处理失败"
            fi
        fi
    done

    shopt -u nullglob

    mark_complete "$checkpoint"
    record_time "步骤8_BiG-Map_map" "完成"
    log_info "步骤8完成: 处理 ${processed} 个样本，跳过 ${skipped} 个样本"
}

# ==========================================
# 步骤9: 合并结果
# ==========================================
step9_merge_results() {
    local step_name="step9_merge_results"
    local checkpoint="$OUTPUT_DIR/0_pipeline_info/checkpoints/${step_name}.complete"

    if is_step_complete "$checkpoint"; then
        log_info "步骤9: 合并结果已完成，跳过"
        return 0
    fi

    record_time "步骤9_合并结果" "开始"
    log_info "开始步骤9: 合并定量结果"

    mkdir -p "$BIGMAP_SUMMARY"

    # 查找所有结果文件
    local result_files=()
    while IFS= read -r -d '' file; do
        result_files+=("$file")
    done < <(find "$BIGMAP_MAP_OUT" -name "BiG-MAP.map.results.ALL.csv" -print0)

    if [[ ${#result_files[@]} -eq 0 ]]; then
        log_warn "未找到定量结果文件"
        return 0
    fi

    # 如果只有一个样本，不需要合并
    if [[ ${#result_files[@]} -eq 1 ]]; then
        log_info "只有一个样本，直接复制结果"
        cp "${result_files[0]}" "$BIGMAP_SUMMARY/merged_BiG-MAP_results_ALL.csv"
        mark_complete "$checkpoint"
        record_time "步骤9_合并结果" "完成"
        log_info "步骤9完成: 单样本结果已复制"
        return 0
    fi

    local output_file="$BIGMAP_SUMMARY/merged_BiG-MAP_results_ALL.csv"

    # 第一个文件完整复制
    cp "${result_files[0]}" "$output_file"

    # 其余文件横向合并
    for f in "${result_files[@]:1}"; do
        paste -d ',' "$output_file" <(cut -d ',' -f 2- "$f") > "${output_file}.tmp"
        mv "${output_file}.tmp" "$output_file"
    done

    mark_complete "$checkpoint"
    record_time "步骤9_合并结果" "完成"
    log_info "步骤9完成: 合并了 ${#result_files[@]} 个样本的结果"
    log_info "合并结果文件: $output_file"
}

# ==========================================
# 主函数
# ==========================================
main() {
    log_info "=========================================="
    log_info "   原核生物BGC分析模块启动"
    log_info "=========================================="

    step1_bakta || exit 1
    step2_antismash || exit 1
    step3_deepbgc || exit 0  # 可选步骤
    step4_extract_deepbgc || exit 0  # 可选步骤
    step5_eggnog_mapper || exit 0  # 可选步骤
    step6_extract_gbk || exit 1
    step7_bigmap_family || exit 1
    step8_bigmap_map || exit 1
    step9_merge_results || exit 0

    log_info "=========================================="
    log_info "   原核生物BGC分析模块完成"
    log_info "=========================================="
}

main "$@"