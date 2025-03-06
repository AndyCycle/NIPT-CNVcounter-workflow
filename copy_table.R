#!/usr/bin/env Rscripts

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
base_input_dir <- "your_input_dir"
hospital_name <- "your_hospital_name"
window_size <- as.numeric(your_window_size)
segs_number <- as.numeric(your_segs_number)

# 验证参数有效性
if (is.na(window_size) || is.na(segs_number)) {
    stop("参数设置错误: window_size 或 segs_number 无效")
}

cat(sprintf("运行参数确认:
"))
cat(sprintf("- 医院: %s
", hospital_name))
cat(sprintf("- 窗口大小: %d
", window_size))
cat(sprintf("- 分段数: %d
", segs_number))

# 根据参数构建目录名称
window_dir_name <- paste0("window_", window_size)
segs_dir <- paste0("segs_", segs_number)

# 输入目录：存放所有 .segments.rds 文件的路径
input_dir <- file.path(base_input_dir, hospital_name, segs_dir, window_dir_name)    #自行修改为符合实际情况的输入目录

# 输出目录：存放最终整合结果（cnv矩阵）的路径
output_base_dir <- "your_output_base_dir"
output_dir <- file.path(output_base_dir, hospital_name, window_dir_name, segs_dir)    #自行修改为符合实际情况的输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义最终输出文件路径（只生成cnv矩阵的RDS文件）
output_rds <- file.path(output_dir, "cnv_table.rds")

# 定义进度跟踪与错误日志文件
progress_file <- file.path(output_dir, "processed_samples.txt")
error_log_file <- file.path(output_dir, "error_log.txt")

# 初始化进度记录文件；若已有内容则读取已处理的样本ID（用于断点续处理）
if (!file.exists(progress_file)) file.create(progress_file)
processed_samples <- if (file.info(progress_file)$size > 0) readLines(progress_file) else character()

# 定义锁文件路径
lock_file <- file.path(output_dir, "cnv_table.lock")

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

# 检查和恢复函数
check_and_recover <- function() {
  if (file.exists(lock_file)) {
    cat("检测到上次处理可能异常中断，进行恢复...\n")
    # 检查临时文件是否存在
    temp_rds <- paste0(output_rds, ".temp")
    temp_progress <- paste0(progress_file, ".temp")

    # 如果临时文件存在，尝试恢复
    if (file.exists(temp_rds)) {
      tryCatch({
        # 尝试读取临时文件验证其完整性
        temp_data <- readRDS(temp_rds)
        # 如果成功读取，用临时文件替换原文件
        file.rename(temp_rds, output_rds)
        cat("成功从临时文件恢复数据\n")
      }, error = function(e) {
        cat("临时文件损坏，将使用原始文件继续处理\n")
        log_error("RECOVERY", sprintf("临时文件恢复失败: %s", e$message))
        if (file.exists(temp_rds)) file.remove(temp_rds)
      })
    }

    # 清理锁文件
    file.remove(lock_file)
  }
}

# 添加时间估算函数
calculate_eta <- function(completed_batches, total_batches, elapsed_time) {
    avg_time_per_batch <- elapsed_time / completed_batches
    remaining_batches <- total_batches - completed_batches
    eta_seconds <- avg_time_per_batch * remaining_batches

    # 转换为可读格式
    hours <- floor(eta_seconds / 3600)
    minutes <- floor((eta_seconds %% 3600) / 60)

    return(sprintf("%dh %dm", hours, minutes))
}

# 主函数（分批次处理）
main <- function() {
  cat("正在生成全局窗口信息...\n")
  region_info <- generate_region_info()

  # 启动处理前检查和恢复
  check_and_recover()

  # 定义批次大小
  batch_size <- 10
  total_files <- length(rds_files)
  batches <- split(rds_files, ceiling(seq_along(rds_files) / batch_size))
  total_batches <- length(batches)
  start_time <- Sys.time()

  # 如果输出文件已存在，检查其中已有的样本
  existing_samples <- character(0)
  if (file.exists(output_rds)) {
    result_df <- readRDS(output_rds)
    existing_samples <- setdiff(names(result_df), c("chr", "start", "end"))
  }

  # 合并两个来源的已处理样本信息
  processed_samples <- unique(c(processed_samples, existing_samples))

  # 处理每个批次
  for (batch_idx in seq_along(batches)) {
    batch <- batches[[batch_idx]]
    current_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(current_time, start_time, units="secs"))

    # 计算预计剩余时间
    if (batch_idx > 1) {
      eta <- calculate_eta(batch_idx - 1, total_batches, elapsed_time)
    } else {
      eta <- "计算中..."
    }

    cat(sprintf("\n========== 批次进度 %d/%d (%.1f%%) ==========\n",
                batch_idx, total_batches, batch_idx / total_batches * 100))
    cat(sprintf("预计剩余时间: %s\n", eta))
    cat(sprintf("当前批次样本数: %d\n", length(batch)))

    # 创建锁文件
    write(Sys.time(), lock_file)

    # 处理批次中的每个样本
    for (rds_file in batch) {
      cat(sprintf("处理文件: %s\n", basename(rds_file)))
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
    }

    # 批次处理完成，强制进行垃圾回收
    cat(sprintf("批次 %d 处理完成，强制进行垃圾回收\n", batch_idx))
    gc()

    # 删除锁文件
    file.remove(lock_file)
  }

  cat(sprintf("所有样本处理完成，最终结果已保存至: %s\n", output_rds))
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

