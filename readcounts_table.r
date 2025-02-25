library(data.table)
library(dplyr)
library(readr)
library(parallel)
library(GenomicRanges)

# 基础参数：输入文件根目录、医院名称、窗口大小和分段文件夹编号
base_input_dir <- "your_base_input_dir"
base_output_dir <- "your_base_output_dir"
hospital_name <- "your_hospital_name"
window_size <- as.numeric(your_window_size)
hmmcopy_number <- as.numeric(your_hmmcopy_number)

# 根据参数构建目录名称
window_dir_name <- paste0("window_", window_size)
hmmcopy_dir <- paste0("hmmcopy_", hmmcopy_number)
readcounts_table_dir <- paste0("readcounts_table_", hmmcopy_number)

# 构建完整的输入输出路径
input_dir <- file.path(base_input_dir, hospital_name, hmmcopy_dir, window_dir_name)
output_dir <- file.path(base_output_dir, hospital_name, window_dir_name, readcounts_table_dir)

# 若输出目录不存在则创建
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义输出文件路径（RDS格式）和进度跟踪文件
output_samples_rds <- file.path(output_dir, "readcounts_table.rds") # 输出的文件名需要根据实际情况修改
progress_file <- file.path(output_dir, "processed_samples.txt") # 样本处理进度跟踪文件

# 函数：读取单个样本文件并提取所需列，同时排除Y染色体
read_sample_data <- function(file_path) {
  tryCatch({
    dt <- readRDS(file_path)
    # 只保留1-22和X染色体的数据
    dt <- dt[grep("^chr(\\d+|X)$", dt$chr), ]
    dt_subset <- dt[, c("chr", "start", "end", "cor.map")]
    return(dt_subset)
  }, error = function(e) {
    message(sprintf("Error reading file %s: %s", file_path, e$message))
    return(NULL)
  })
}

# 辅助函数：对染色体进行自然排序，返回排序后的索引
chr_sort <- function(chr_vec) {
  valid_idx <- grep("^chr(\\d+|X)$", chr_vec)
  chr_valid <- chr_vec[valid_idx]
  nums <- gsub("chr", "", chr_valid)
  nums[nums == "X"] <- "23"
  nums <- as.numeric(nums)
  idx_order <- order(nums)
  return(valid_idx[idx_order])
}

# 修改创建初始 RDS 文件函数
create_initial_rds <- function(first_sample_dt, output_file) {
  sort_idx <- chr_sort(first_sample_dt$chr)
  first_sample_dt <- first_sample_dt[sort_idx, ]

  n_windows <- nrow(first_sample_dt)
  sample_id <- sub("_[0-9]+bp\\.correctedReadcount\\.rds$", "", basename(rds_files[1]))
  readcounts <- first_sample_dt$cor.map

  dt_table <- data.table(
    chr = first_sample_dt$chr,
    start = first_sample_dt$start,
    end = first_sample_dt$end
  )

  # 直接使用样本ID作为列名
  dt_table[[sample_id]] <- readcounts

  saveRDS(dt_table, output_file)
  return(n_windows)
}

# 添加临时文件和锁文件机制
update_rds_with_sample <- function(rds_file, sample_vector, sample_id) {
  # 定义临时文件和锁文件
  temp_rds <- paste0(rds_file, ".temp")
  lock_file <- paste0(rds_file, ".lock")

  tryCatch({
    # 创建锁文件
    write(Sys.time(), lock_file)

    # 读取原始数据
    dt_table <- readRDS(rds_file)

    if (length(sample_vector) != nrow(dt_table)) {
      message(sprintf("样本 %s 行数与窗口信息不匹配，跳过该样本。", sample_id))
      return(FALSE)
    }

    # 更新数据
    dt_table[[sample_id]] <- sample_vector

    # 先写入临时文件
    saveRDS(dt_table, temp_rds)

    # 使用原子操作替换原文件
    file.rename(temp_rds, rds_file)

    # 清理内存
    rm(dt_table)
    gc()

    return(TRUE)
  }, error = function(e) {
    message(sprintf("处理样本 %s 时发生错误: %s", sample_id, e$message))
    return(FALSE)
  }, finally = {
    # 清理临时文件和锁文件
    if (file.exists(temp_rds)) file.remove(temp_rds)
    if (file.exists(lock_file)) file.remove(lock_file)
  })
}

# 检查和恢复函数
check_and_recover <- function(rds_file) {
  lock_file <- paste0(rds_file, ".lock")
  temp_rds <- paste0(rds_file, ".temp")

  # 检查是否存在锁文件，表明上次处理可能异常中断
  if (file.exists(lock_file)) {
    message("检测到上次处理可能异常中断，进行恢复...")

    # 如果存在临时文件，检查其完整性
    if (file.exists(temp_rds)) {
      tryCatch({
        # 尝试读取临时文件验证其完整性
        temp_data <- readRDS(temp_rds)
        # 如果成功读取，用临时文件替换原文件
        file.rename(temp_rds, rds_file)
        message("成功从临时文件恢复数据")
      }, error = function(e) {
        message("临时文件损坏，将使用原始文件继续处理")
        if (file.exists(temp_rds)) file.remove(temp_rds)
      })
    }

    file.remove(lock_file)
  }
}

# 修改批处理函数
process_samples_in_batches <- function(rds_files, batch_size = 10) {
  total_files <- length(rds_files)
  batches <- split(rds_files, ceiling(seq_along(rds_files)/batch_size))

  # 启动处理前检查和恢复
  check_and_recover(output_samples_rds)

  for (batch in batches) {
    # 处理每个批次的样本
    for (file in batch) {
      sample_id <- sub("_[0-9]+bp\\.correctedReadcount\\.rds$", "", basename(file))

      # 检查是否已处理
      if (sample_id %in% processed_samples) {
        message(sprintf("跳过已处理样本: %s", sample_id))
        next
      }

      message(sprintf("开始处理样本: %s", sample_id))

      sample_dt <- read_sample_data(file)
      if (is.null(sample_dt)) {
        message(sprintf("跳过处理失败的样本: %s", sample_id))
        next
      }

      if (nrow(sample_dt) != n_windows) {
        message(sprintf("样本 %s 行数与窗口信息不匹配，跳过该样本。", sample_id))
        next
      }

      # 更新数据并记录进度
      if (update_rds_with_sample(output_samples_rds, sample_dt$cor.map, sample_id)) {
        write(sample_id, file = progress_file, append = TRUE)
        processed_samples <- c(processed_samples, sample_id)
        message(sprintf("已处理样本 %s (%d/%d)", sample_id, length(processed_samples), total_files))
      }

      rm(sample_dt)
      gc()
    }
    gc()
  }
}

# 列出所有符合条件的 RDS 文件
rds_files <- list.files(path = input_dir,
                        pattern = "\\.correctedReadcount\\.rds$",
                        full.names = TRUE)
if (length(rds_files) == 0) {
  stop("在输入目录中未找到任何 RDS 文件")
}

# 加载已处理样本记录（若存在进度文件）
if (file.exists(progress_file)) {
  processed_samples <- readLines(progress_file)
} else {
  processed_samples <- character(0)
}

total_files <- length(rds_files)

# 如果 RDS 输出文件不存在，则用第一个样本数据构建初始数据表
if (!file.exists(output_samples_rds)) {
  message("处理第一个文件以获取窗口信息及初始 readcounts 数据...")
  first_sample_dt <- read_sample_data(rds_files[1])
  if (is.null(first_sample_dt)) {
    stop("第一个样本文件处理失败")
  }
  n_windows <- create_initial_rds(first_sample_dt, output_samples_rds)

  sample_id_first <- sub("_[0-9]+bp\\.correctedReadcount\\.rds$", "", basename(rds_files[1]))
  write(sample_id_first, file = progress_file)
  processed_samples <- c(processed_samples, sample_id_first)
  message(sprintf("已处理样本 %s (1/%d)", sample_id_first, total_files))
} else {
  dt_table <- readRDS(output_samples_rds)
  n_windows <- nrow(dt_table)
  rm(dt_table) # 立即释放内存
  gc()
}

# 使用批处理函数处理样本
process_samples_in_batches(rds_files, batch_size = 10)

message("所有样本处理完成！")
message(sprintf("共处理样本数量: %d", length(processed_samples)))
message("输出文件: ", output_samples_rds)
message("进度跟踪文件: ", progress_file)

