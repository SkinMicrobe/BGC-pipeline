#!/bin/bash

set -euo pipefail

# ==========================================
# BGC Pipeline 主控脚本
# ==========================================

# 脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$SKILL_ROOT/config"

# 加载默认参数
DEFAULT_PARAMS="$CONFIG_DIR/default_params.yaml"
USER_INPUT="$CONFIG_DIR/user_input.yaml"

# 日志函数
log_info() { echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_warn() { echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2; }
log_error() { echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') $*" >&2; }

# ==========================================
# 解析YAML配置（简化版）
# ==========================================
parse_yaml() {
    local file="$1"
    local prefix="$2"
    grep -E "^[[:space:]]*[^#]" "$file" | \
    sed "s/^\([a-z_]*\)[[:space:]]*:[[:space:]]*\(.*\)/\1=\2/" | \
    while read -r line; do
        key="${line%%=*}"
        val="${line#*=}"
        # 去除引号
        val="${val%\"}"
        val="${val#\"}"
        echo "${prefix}${key}=${val}"
    done
}

# ==========================================
# 交互式输入函数
# ==========================================
prompt_user_input() {
    echo "=========================================="
    echo "   BGC Pipeline - 生物合成基因簇分析流程"
    echo "=========================================="
    echo

    # 模式选择
    while true; do
        read -p "请选择模式 [1.原核生物 prokaryote, 2.真核生物 eukaryote]: " mode_choice
        case $mode_choice in
            1) MODE="prokaryote"; break ;;
            2) MODE="eukaryote"; break ;;
            prokaryote) MODE="prokaryote"; break ;;
            eukaryote) MODE="eukaryote"; break ;;
            *) echo "无效选择，请重新输入" ;;
        esac
    done

    # 输入目录
    while true; do
        read -p "请输入FASTA/MAG文件目录: " INPUT_DIR
        if [[ -d "$INPUT_DIR" ]]; then
            break
        else
            echo "目录不存在，请重新输入"
        fi
    done

    # 输出目录
    read -p "请输入输出目录（将自动创建）: " OUTPUT_DIR

    # 线程数
    read -p "请输入线程数 [默认32]: " THREADS_INPUT
    THREADS=${THREADS_INPUT:-32}

    # 原核模式额外参数
    if [[ "$MODE" == "prokaryote" ]]; then
        while true; do
            read -p "请输入测序fastq文件目录（用于定量）: " FASTQ_DIR
            if [[ -d "$FASTQ_DIR" ]]; then
                break
            else
                echo "目录不存在，请重新输入"
            fi
        done

        read -p "是否运行deepBGC? [y/N]: " RUN_DEEPBGC
        read -p "是否运行eggNOG? [y/N]: " RUN_EGGNOG
    fi

    # 可选步骤（原核）
    read -p "是否跳过已完成的步骤（断点续传）? [Y/n]: " SKIP_EXISTING
    SKIP_EXISTING=${SKIP_EXISTING:-y}
}

# ==========================================
# 初始化Conda
# ==========================================
init_conda() {
    log_info "初始化Conda环境..."
    if ! command -v conda &> /dev/null; then
        log_error "未找到conda命令"
        exit 1
    fi
    eval "$(conda shell.bash hook)"
}

# ==========================================
# 激活环境检查
# ==========================================
check_conda_env() {
    local env_name="$1"
    log_info "检查Conda环境: $env_name"
    if conda env list | grep -q "^${env_name} "; then
        log_info "环境 $env_name 存在"
        return 0
    else
        log_warn "环境 $env_name 不存在，请先安装"
        return 1
    fi
}

# ==========================================
# 创建输出目录结构
# ==========================================
create_output_structure() {
    log_info "创建输出目录结构..."

    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR/0_pipeline_info/checkpoints"
    mkdir -p "$OUTPUT_DIR/0_pipeline_info/logs"

    if [[ "$MODE" == "prokaryote" ]]; then
        mkdir -p "$OUTPUT_DIR/prok_results/1_bakta_annotation"
        mkdir -p "$OUTPUT_DIR/prok_results/2_antismash_out"
        mkdir -p "$OUTPUT_DIR/prok_results/3_deepbgc_out"
        mkdir -p "$OUTPUT_DIR/prok_results/4_deepbgc_extracted"
        mkdir -p "$OUTPUT_DIR/prok_results/5_eggnog_mapper_out"
        mkdir -p "$OUTPUT_DIR/prok_results/6_bigmap_family/aggregated_gbks"
        mkdir -p "$OUTPUT_DIR/prok_results/6_bigmap_family/bigscape_mash_results"
        mkdir -p "$OUTPUT_DIR/prok_results/7_bigmap_map_quant"
        mkdir -p "$OUTPUT_DIR/prok_results/8_bigmap_summary"
    else
        mkdir -p "$OUTPUT_DIR/euk_results/1_metaeuk_annotation"
        mkdir -p "$OUTPUT_DIR/euk_results/2_antismash_out"
        mkdir -p "$OUTPUT_DIR/euk_results/3_gbk_extracted"
        mkdir -p "$OUTPUT_DIR/euk_results/4_bigmap_family"
    fi

    # 创建临时目录
    mkdir -p "$TMP_DIR"
}

# ==========================================
# 保存运行配置
# ==========================================
save_config() {
    local config_file="$OUTPUT_DIR/0_pipeline_info/run_config.log"
    log_info "保存运行配置到 $config_file"

    cat > "$config_file" << EOF
# BGC Pipeline 运行配置
# 生成时间: $(date '+%Y-%m-%d %H:%M:%S')

# ============================================
# 基本信息
# ============================================
MODE=$MODE
INPUT_DIR=$INPUT_DIR
OUTPUT_DIR=$OUTPUT_DIR
THREADS=$THREADS
TMP_DIR=$TMP_DIR

# ============================================
# 原核模式参数
# ============================================
FASTQ_DIR=${FASTQ_DIR:-}
RUN_DEEPBGC=$([[ "${RUN_DEEPBGC:-n}" =~ ^[Yy] ]] && echo "true" || echo "false")
RUN_EGGNOG=$([[ "${RUN_EGGNOG:-n}" =~ ^[Yy] ]] && echo "true" || echo "false")

# ============================================
# 数据库路径
# ============================================
BAKTA_DB=$BAKTA_DB
PFAM_HMM=$PFAM_HMM
EGGNOG_DB=$EGGNOG_DB
UNIREF50_DB=$UNIREF50_DB

# ============================================
# 断点续传设置
# ============================================
SKIP_EXISTING=$([[ "${SKIP_EXISTING:-y}" =~ ^[Yy] ]] && echo "true" || echo "false")
EOF
}

# ==========================================
# 检查数据库路径
# ==========================================
check_databases() {
    log_info "检查数据库路径..."

    local check_failed=false

    if [[ "$MODE" == "prokaryote" ]]; then
        if [[ ! -d "$BAKTA_DB" ]]; then
            log_error "Bakta数据库不存在: $BAKTA_DB"
            check_failed=true
        fi
        if [[ ! -f "$PFAM_HMM" ]]; then
            log_error "Pfam数据库不存在: $PFAM_HMM"
            check_failed=true
        fi
    else
        if [[ ! -d "$UNIREF50_DB" ]]; then
            log_error "UniRef50数据库不存在: $UNIREF50_DB"
            check_failed=true
        fi
    fi

    if [[ "$check_failed" == "true" ]]; then
        log_error "数据库检查失败，请检查config/default_params.yaml中的路径"
        exit 1
    fi

    log_info "数据库检查通过"
}

# ==========================================
# 运行时间记录
# ==========================================
record_time() {
    local step="$1"
    local action="$2"
    local time_file="$OUTPUT_DIR/0_pipeline_info/run_time.log"

    if [[ "$action" == "start" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $step - 开始" >> "$time_file"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $step - 完成" >> "$time_file"
    fi
}

# ==========================================
# 检查步骤是否完成
# ==========================================
is_step_complete() {
    local marker="$1"
    if [[ "${SKIP_EXISTING:-true}" == "true" ]] && [[ -f "$marker" ]]; then
        return 0
    fi
    return 1
}

# ==========================================
# 主流程
# ==========================================
main() {
    echo
    echo "=========================================="
    echo "   BGC Pipeline - 生物合成基因簇分析流程"
    echo "=========================================="
    echo

    # 检查是否已有配置文件
    if [[ -f "$USER_INPUT" ]] && grep -q "mode:" "$USER_INPUT"; then
        log_info "检测到已有配置文件: $USER_INPUT"
        read -p "是否使用现有配置？[Y/n]: " USE_EXISTING
        USE_EXISTING=${USE_EXISTING:-y}
    fi

    if [[ "${USE_EXISTING:-n}" =~ ^[Yy] ]]; then
        # 从配置文件读取（简化版，实际应用需要更完整的YAML解析）
        log_info "从配置文件读取参数..."
        source "$USER_INPUT" 2>/dev/null || true
        # 如果配置文件解析失败，重新输入
        if [[ -z "${MODE:-}" ]]; then
            prompt_user_input
        fi
    else
        prompt_user_input
    fi

    # 读取默认参数中的数据库路径
    while IFS='=' read -r key val; do
        case "$key" in
            bakta_db) BAKTA_DB="$val" ;;
            pfam_hmm) PFAM_HMM="$val" ;;
            eggnog_db) EGGNOG_DB="$val" ;;
            uniref50_db) UNIREF50_DB="$val" ;;
            tmp_dir) TMP_DIR="${val:-/tmp/bgcpipeline_$$}" ;;
        esac
    done < <(parse_yaml "$DEFAULT_PARAMS" "")

    # 初始化
    init_conda

    # 检查Conda环境
    log_info "检查所需Conda环境..."
    local envs_needed="bakta antismash bigmap"
    if [[ "$MODE" == "eukaryote" ]]; then
        envs_needed="$envs_needed meta"
    fi
    # 原核模式可选环境
    if [[ "$MODE" == "prokaryote" ]]; then
        if [[ "${RUN_DEEPBGC:-n}" =~ ^[Yy] ]]; then
            envs_needed="$envs_needed deepbgc"
        fi
        if [[ "${RUN_EGGNOG:-n}" =~ ^[Yy] ]]; then
            envs_needed="$envs_needed emapper"
        fi
    fi

    for env in $envs_needed; do
        check_conda_env "$env" || exit 1
    done

    # 检查数据库
    check_databases

    # 创建输出目录
    create_output_structure

    # 保存配置
    save_config

    echo
    log_info "配置完成，即将开始运行..."
    log_info "模式: $MODE"
    log_info "输入目录: $INPUT_DIR"
    log_info "输出目录: $OUTPUT_DIR"
    log_info "线程数: $THREADS"
    echo

    # 运行对应模式的模块
    if [[ "$MODE" == "prokaryote" ]]; then
        log_info "开始运行原核生物模块..."
        bash "$SCRIPT_DIR/module_prok.sh" || exit 1
    else
        log_info "开始运行真核生物模块..."
        bash "$SCRIPT_DIR/module_euk.sh" || exit 1
    fi

    echo
    log_info "=========================================="
    log_info "   BGC Pipeline 运行完成！"
    log_info "=========================================="
    log_info "结果目录: $OUTPUT_DIR"
    log_info "日志目录: $OUTPUT_DIR/0_pipeline_info/logs"
    echo
}

# 运行主函数
main "$@"