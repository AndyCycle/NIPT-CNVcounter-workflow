#!/usr/bin/env Rscript

# 设置工作目录（根据实际情况修改）
setwd("your_working_directory")

# 定义参数：可以根据需要调整医院名称、窗口大小和分段数
hospital_names <- c("your_hospital_name") # 医院名称
window_sizes <- c(your_window_size) # 窗口大小
segs_count <- 12 # 文件夹分段数

# 模板脚本内容（新版本处理脚本，无 output_regions.csv，同时实时输出进度信息）
template_script <- paste0(
'#!/usr/bin/env Rscripts

# ---------------------------------------------------------------------------
# 加载必要的库
library(data.table)
library(GenomicRanges)
library(dplyr)
library(readr)
library(gtools)  # 用于染色体自然排序

# ---------------------------------------------------------------------------
# 1. 设置基础参数 & 构建输入/输出目录

# 基础参数：输入文件根目录、医院名称、窗口大小和分段文件夹编号
base_input_dir <- "your_base_input_dir"
hospital_name <- {{hospital_name}}
window_size <- as.numeric({{window_size}})
segs_number <- as.numeric({{segs_number}})

# 验证参数有效性
if (is.na(window_size) || is.na(segs_number)) {
    stop("参数设置错误: window_size 或 segs_number 无效")
}

cat(sprintf("运行参数确认:\n"))
cat(sprintf("- 医院: %s\n", hospital_name))
cat(sprintf("- 窗口大小: %d\n", window_size))
cat(sprintf("- 分段数: %d\n", segs_number))

# 根据参数构建目录名称
window_dir_name <- paste0("window_", window_size)
segs_dir <- paste0("segs_", segs_number)

# 输入目录：存放所有 .segments.rds 文件的路径
input_dir <- file.path(base_input_dir, hospital_name, segs_dir, window_dir_name)

# 输出目录：存放最终整合结果（cnv矩阵）的路径
output_base_dir <- "your_output_base_dir"
output_dir <- file.path(output_base_dir, hospital_name, window_dir_name, segs_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义最终输出文件路径（只生成cnv矩阵的RDS文件）
output_rds <- file.path(output_dir, "cnv_table.rds")

# 定义进度跟踪与错误日志文件
progress_file <- file.path(output_dir, "processed_samples.txt")
error_log_file <- file.path(output_dir, "error_log.txt")

# 初始化进度记录文件；若已有内容则读取已处理的样本ID（用于断点续处理）
if (!file.exists(progress_file)) file.create(progress_file)
processed_samples <- if (file.info(progress_file)$size > 0) readLines(progress_file) else character()

# ---------------------------------------------------------------------------
# 获取所有符合条件的 .segments.rds 文件
rds_files <- list.files(input_dir, pattern = "[.]segments[.]rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("未找到任何 .segments.rds 文件。")

# ---------------------------------------------------------------------------
# 2. 定义辅助函数

# 日志记录函数
log_error <- function(file, message, append_timestamp = TRUE) {
    error_msg <- if (append_timestamp) {
        sprintf("[%s] %s: %s
", format(Sys.time()), basename(file), message)
    } else {
        sprintf("%s: %s
", basename(file), message)
    }
    cat(error_msg, file = error_log_file, append = TRUE)
}

# 提取样本信息函数（解析文件名获得 sample_id 和窗口大小信息）
extract_sample_info <- function(file_path) {
    tryCatch({
        name_parts <- strsplit(basename(file_path), "[.]")[[1]][1]
        parts <- strsplit(name_parts, "_")[[1]]
        # 获取窗口大小部分，并去掉 "bp" 后缀
        window_size_str <- gsub("bp$", "", parts[length(parts)])
        window_size_val <- as.numeric(window_size_str)
        # 拼接其余部分作为 sample_id
        sample_id <- paste(parts[-length(parts)], collapse = "_")
        list(
            sample_id = sample_id,
            window_size = window_size_val,
            file_path = file_path
        )
    }, error = function(e) {
        log_error(file_path, "样本信息提取失败")
        return(NULL)
    })
}

# 生成全局窗口信息函数
# 遍历所有样本，获取每条染色体的最大结束位置，并根据 window_size 生成统一区域信息（chr, start, end）
generate_region_info <- function() {
    max_chr_lengths <- list()

    for (file in rds_files) {
        sample_info <- extract_sample_info(file)
        if (is.null(sample_info)) next
        if (sample_info$window_size != window_size) next

        tryCatch({
            rds_data <- readRDS(file)
            # 合并多个染色体的数据，添加新列 "chr" 标识染色体名
            combined_data <- rbindlist(rds_data, idcol = "chr")
            combined_data[, chr := as.character(chr)]

            for (chr_name in unique(combined_data$chr)) {
                current_max <- max(combined_data[chr == chr_name, end], na.rm = TRUE)
                if (is.null(max_chr_lengths[[chr_name]]) ||
                    current_max > max_chr_lengths[[chr_name]]) {
                    max_chr_lengths[[chr_name]] <- current_max
                }
            }
        }, error = function(e) {
            log_error(file, e$message)
            next
        })
    }

    if (length(max_chr_lengths) == 0) stop("无法收集染色体长度信息。")

    # 自然排序染色体
    sorted_chr <- mixedsort(names(max_chr_lengths))

    regions_list <- lapply(sorted_chr, function(chr) {
        chr_length <- max_chr_lengths[[chr]]
        num_windows <- ceiling(chr_length / window_size)
        data.frame(
            chr = chr,
            start = ((0:(num_windows - 1)) * window_size) + 1,
            end = pmin((1:num_windows) * window_size, chr_length),
            stringsAsFactors = FALSE
        )
    })

    region_info <- do.call(rbind, regions_list)
    return(region_info)
}

# 根据单个样本的 segmentation 数据将拷贝数映射到全局窗口上
process_sample_data <- function(rds_file, region_info) {
    tryCatch({
        sample_info <- extract_sample_info(rds_file)
        if (is.null(sample_info)) return(NULL)
        sample_id <- sample_info$sample_id

        # 如果该样本已处理，跳过（用于断点续处理）
        if (sample_id %in% processed_samples) {
            cat(sprintf("样本 %s 已处理，跳过
", sample_id))
            return(NULL)
        }

        # 读取样本分段数据，已预处理为固定窗口倍数，列表中每个元素对应一条染色体
        seg_data <- readRDS(rds_file)
        # 初始化长度与全局窗口数相同的向量保存拷贝数数据
        sample_cnv <- rep(NA, nrow(region_info))

        # 逐个染色体读取和处理
        seg_data <- readRDS(rds_file)
        for (chr in names(seg_data)) {
            chr_seg <- seg_data[[chr]]
            chr_windows_idx <- which(region_info$chr == chr)
            if (length(chr_windows_idx) == 0) next

            chr_windows <- region_info[chr_windows_idx, ]

            # 使用数据表进行更高效的操作
            chr_seg_dt <- as.data.table(chr_seg)
            chr_windows_dt <- as.data.table(chr_windows)

            # 对每个窗口，找到对应的片段
            for (i in seq_len(nrow(chr_windows_dt))) {
                window_start <- chr_windows_dt$start[i]
                matching_seg <- chr_seg_dt[window_start >= start & window_start <= end]
                if (nrow(matching_seg) > 0) {
                    cnv_value <- as.numeric(matching_seg$state[1]) - 1
                    sample_cnv[chr_windows_idx[i]] <- cnv_value
                }
            }

            # 及时清理临时对象
            rm(chr_seg, chr_windows, chr_seg_dt, chr_windows_dt)
            gc()
        }
        rm(seg_data)
        gc()

        return(list(
            sample_id = sample_id,
            cnv_data = sample_cnv
        ))

    }, error = function(e) {
        log_error(rds_file, e$message)
        return(NULL)
    })
}

# ---------------------------------------------------------------------------
# 3. 主函数（边处理边写入）
#    - 1) 首先生成全局窗口信息 region_info
#    - 2) 对每个样本调用 process_sample_data 得到该样本的拷贝数据向量
#    - 3) 每处理完一个样本，就从磁盘读入现有 rds 文件（若不存在则用 region_info 初始化），
#         追加当前样本的列后更新写回，形成边处理边写入的机制
main <- function() {
    cat("正在生成全局窗口信息...
")
    region_info <- generate_region_info()

    cat("开始边处理边写入样本数据...
")
    chunk_size <- 10 # 每处理10个样本强制进行一次垃圾回收

    # 如果输出文件已存在，检查其中已有的样本
    existing_samples <- character(0)
    if (file.exists(output_rds)) {
        result_df <- readRDS(output_rds)
        existing_samples <- setdiff(names(result_df), c("chr", "start", "end"))
    }

    # 合并两个来源的已处理样本信息
    processed_samples <- unique(c(processed_samples, existing_samples))

    for (i in seq_along(rds_files)) {
        rds_file <- rds_files[i]
        cat(sprintf("处理文件: %s
", basename(rds_file)))
        sample_result <- process_sample_data(rds_file, region_info)

        if (!is.null(sample_result)) {
            # 添加事务性写入
            temp_rds <- paste0(output_rds, ".temp")
            temp_progress <- paste0(progress_file, ".temp")

            tryCatch({
                # 保存数据到临时文件
                result_df <- if (file.exists(output_rds)) readRDS(output_rds) else region_info
                result_df[[sample_result$sample_id]] <- sample_result$cnv_data
                saveRDS(result_df, temp_rds)
                write(sample_result$sample_id, temp_progress)

                # 如果成功，替换原文件
                file.rename(temp_rds, output_rds)
                file.append(progress_file, temp_progress)
                file.remove(temp_progress)
            }, error = function(e) {
                # 清理临时文件
                if (file.exists(temp_rds)) file.remove(temp_rds)
                if (file.exists(temp_progress)) file.remove(temp_progress)
                stop(sprintf("保存样本 %s 时发生错误: %s", sample_result$sample_id, e$message))
            })
        }

        # 每处理chunk_size个样本强制进行一次垃圾回收
        if (i %% chunk_size == 0) {
            gc()
        }
    }

    cat(sprintf("所有样本处理完成，最终结果已保存至: %s
", output_rds))
}

# ---------------------------------------------------------------------------
# 4. 程序执行入口
if (!interactive()) {
    tryCatch({
        main()
    }, error = function(e) {
        log_error("MAIN", sprintf("程序执行失败: %s", e$message))
        quit(status = 1)
    })
}
'
, collapse = "\n")

for (hospital in hospital_names) {
    for (window in window_sizes) {
        for (seg in 1:segs_count) {
            # 替换模板中的占位符
            script_content <- template_script
            script_content <- gsub('{{hospital_name}}', sprintf('"%s"', hospital), script_content, fixed = TRUE)
            script_content <- gsub('{{window_size}}', window, script_content, fixed = TRUE)
            script_content <- gsub('{{segs_number}}', seg, script_content, fixed = TRUE)

            # 构造脚本存放目录
            script_dir <- file.path("sbatch_jobs/scr/segs_table_new", hospital, #脚本存放目录 （根据实际文件夹存放情况修改）
                                    paste0("window_", window),
                                    paste0("segs_", seg))
            dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

            # 定义生成的处理脚本文件名
            script_filename <- file.path(script_dir, "process_segs.R")
            writeLines(script_content, script_filename)
            system(paste("chmod +x", script_filename))

            # 生成对应的 sbatch 脚本
            sbatch_content <- sprintf('#!/bin/bash
#SBATCH --job-name=%s_w%d_s%d
#SBATCH --output=%s_w%d_s%d_slurm_%%j.out
#SBATCH --error=%s_w%d_s%d_slurm_%%j.err
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=150:00:00
#SBATCH --mem=8G

echo "cnv_table process started at: $(date)"

export PATH="your_conda_path/bin:$PATH"  #注意更改为实际conda路径

source activate hmmcopy  # 激活你的conda环境（根据实际情况修改）

# 使用完整的绝对路径
Rscript your_working_directory/sbatch_jobs/scr/segs_table_new/%s/window_%d/segs_%d/process_segs.R

echo "Process completed at: $(date)"
', hospital, window, seg, hospital, window, seg, hospital, window, seg, hospital, window, seg)

            sbatch_filename <- file.path(script_dir, "run_process.sh")
            writeLines(sbatch_content, sbatch_filename)
            system(paste("chmod +x", sbatch_filename))

            cat(sprintf("Generated scripts for %s, window %d, segs %d\n", hospital, window, seg))
        }
    }
}

# 创建一个批量提交所有作业的脚本（与之前类似）
submit_all_content <- '#!/bin/bash

BASE_DIR="your_working_directory" # 工作目录

for hospital in your_hospital_name; do  # 根据实际情况修改hospital_name——你的医院名称（可多个）
    for window in your_window_size; do  # 根据实际情况修改window_size——你的窗口大小
        for seg in $(seq 1 your_segs_count); do  # 根据实际情况修改segs_count——你的文件夹分段数
            script_path="${BASE_DIR}/sbatch_jobs/scr/segs_table_new/${hospital}/window_${window}/segs_${seg}/run_process.sh"
            if [ -f "$script_path" ] && [ -x "$script_path" ]; then
                sbatch ${script_path}
                echo "成功提交作业: ${script_path}"
                sleep 1
            else
                echo "错误: 文件不存在或没有执行权限: ${script_path}"
            fi
        done
    done
done
'

submit_all_filename <- "sbatch_jobs/scr/segs_table_new/submit_all_jobs.sh"
writeLines(submit_all_content, submit_all_filename)
system(paste("chmod +x", submit_all_filename))

cat("Script generation completed!\n")
cat("To submit all jobs, run: ./submit_all_jobs.sh\n")