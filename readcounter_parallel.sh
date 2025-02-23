#!/bin/bash
#SBATCH --job-name=ReadCounts_Parallel
#SBATCH --output=%j_ReadCounts_Parallel.out
#SBATCH --error=%j_ReadCounts_Parallel.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=80:00:00

export PATH="your_conda_path/bin:$PATH"

# 创建临时文件来跟踪处理状态
TEMP_DIR=$(mktemp -d)
FAILED_TASKS="${TEMP_DIR}/failed_tasks.txt"
COMPLETED_TASKS="${TEMP_DIR}/completed_tasks.txt"
touch "${FAILED_TASKS}" "${COMPLETED_TASKS}"

# 激活环境
source activate hmmcopy
if [ $? -ne 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Failed to activate conda environment"  >&2
    exit 1
fi

# 检查 parallel
if ! command -v parallel &> /dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] GNU Parallel not found. Installing via conda..."
    conda install -y -c conda-forge parallel
    if [ $? -ne 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Failed to install GNU Parallel"  >&2
        exit 1
    fi
fi

# 设置变量
READ_COUNTER="your_Read_Counter_path" # hmmcopy_utils/bin/readCounter工具
BAM_DIR="your_bam_dir"
OUTPUT_DIR="your_output_dir"
WINDOW_SIZES=(your_window_size)

mkdir -p "${OUTPUT_DIR}"

# 修改获取BAM文件的方式，使用find命令递归查找所有子目录中的BAM文件
mapfile -t SAMPLE_BAMS < <(find -L "${BAM_DIR}" -type f -name "*.bam")
TOTAL_SAMPLES=${#SAMPLE_BAMS[@]}

if [ ${TOTAL_SAMPLES} -eq 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: No BAM files found in ${BAM_DIR}"  >&2
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found ${TOTAL_SAMPLES} BAM files"

# 改进的处理函数，修改样本ID提取方式
process_readcounts() {
    local BAM_FILE=$1
    local WINDOW=$2
    local RETRY_COUNT=3
    local retry=0

    # 从BAM文件路径中提取样本ID
    SAMPLE_ID=$(basename "${BAM_FILE}" | sed 's/\.sorted\.rmdup\.realign\.BQSR\.bam$//' | cut -d. -f1)

    WINDOW_DIR="${OUTPUT_DIR}/window_${WINDOW}"
    mkdir -p "${WINDOW_DIR}"

    READCOUNTER_OUTPUT="${WINDOW_DIR}/${SAMPLE_ID}_${WINDOW}bp.readcounts.wig"
    TASK_ID="${SAMPLE_ID}_${WINDOW}"

    # 检查是否已经成功处理过
    if [ -f "${READCOUNTER_OUTPUT}" ] && [ -s "${READCOUNTER_OUTPUT}" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] File ${READCOUNTER_OUTPUT} already exists and is non-empty, skipping..."
        echo "${TASK_ID}" >> "${COMPLETED_TASKS}"
        return 0
    fi

    while [ $retry -lt $RETRY_COUNT ]; do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${SAMPLE_ID} with window size ${WINDOW} (attempt $((retry+1))/${RETRY_COUNT})"

        ${READ_COUNTER} -w "${WINDOW}" -q 20 "${BAM_FILE}" | sed '/^fixedStep chrom=chrM/,$d' > "${READCOUNTER_OUTPUT}.tmp"  # -q 20 表示将质量值低于20的读段过滤掉  # sed '/^fixedStep chrom=chrM/,$d' 表示删除chrM染色体后的内容

        if [ $? -eq 0 ] && [ -s "${READCOUNTER_OUTPUT}.tmp" ]; then
            mv "${READCOUNTER_OUTPUT}.tmp" "${READCOUNTER_OUTPUT}"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Successfully processed ${TASK_ID}"
            echo "${TASK_ID}" >> "${COMPLETED_TASKS}"
            return 0
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Attempt $((retry+1)) failed for ${TASK_ID}"  >&2
            rm -f "${READCOUNTER_OUTPUT}.tmp"
            retry=$((retry+1))
            sleep 1
        fi
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Failed to process ${TASK_ID} after ${RETRY_COUNT} attempts"  >&2
    echo "${TASK_ID}" >> "${FAILED_TASKS}"
    return 1
}

export -f process_readcounts
export READ_COUNTER OUTPUT_DIR COMPLETED_TASKS FAILED_TASKS

# 使用parallel执行任务，并设置超时 (自行修改timeout时间和并行jobs数量)
parallel --timeout 3600 --jobs 5 process_readcounts {} ::: "${SAMPLE_BAMS[@]}" ::: "${WINDOW_SIZES[@]}"  # --timeout 3600 表示设置超时时间为3600秒  # --jobs 5 表示同时运行5个任务

# 验证处理结果
TOTAL_EXPECTED=$((TOTAL_SAMPLES * ${#WINDOW_SIZES[@]}))
COMPLETED_COUNT=$(wc -l < "${COMPLETED_TASKS}")
FAILED_COUNT=$(wc -l < "${FAILED_TASKS}")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing summary:"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expected tasks: ${TOTAL_EXPECTED}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed tasks: ${COMPLETED_COUNT}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed tasks: ${FAILED_COUNT}"

if [ -s "${FAILED_TASKS}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed tasks list:"
    cat "${FAILED_TASKS}" | while read -r task; do
        echo "  - ${task}"
    done
fi

# 检查是否所有任务都已完成
if [ "$((COMPLETED_COUNT + FAILED_COUNT))" -ne "${TOTAL_EXPECTED}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Some tasks were not processed. Please check the logs."  >&2
    exit 1
fi

# 清理临时文件
rm -rf "${TEMP_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Process completed"