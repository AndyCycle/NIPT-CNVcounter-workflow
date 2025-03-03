library(HMMcopy)
library(data.table)
library(parallel)

# 设置输入输出目录
input_dir <- "your_input_dir"
output_dir <- "your_output_dir"
progress_dir <- file.path(output_dir, "progress")  # 将progress目录放在output_dir下

# 创建输出目录
dirs <- c(output_dir, progress_dir)
for(dir in dirs) {
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)  # 添加recursive参数确保父目录存在
    }
}

# 从输入目录获取所有待处理样本
get_all_samples <- function() {
    input_files <- list.files(input_dir, pattern="*_10000bp.correctedReadcount.rds$")
    return(gsub("_10000bp.correctedReadcount.rds$", "", input_files))
}

# 检查已处理的样本
get_processed_samples <- function() {
    processed_files <- list.files(output_dir, pattern="*_10000bp.segments.rds$")
    return(gsub("_10000bp.segments.rds$", "", processed_files))
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

# 获取待处理样本
all_samples <- get_all_samples()
processed_samples <- get_processed_samples()
remaining_samples <- setdiff(all_samples, processed_samples)

print(sprintf("总样本数: %d", length(all_samples)))
print(sprintf("已处理样本数: %d", length(processed_samples)))
print(sprintf("待处理样本数: %d", length(remaining_samples)))

# 定义处理单个样本的函数
process_sample <- function(sample_id) {
    print(sprintf("[%s] 开始处理样本: %s",
                 format(Sys.time(), "%H:%M:%S"),
                 sample_id))

    tryCatch({
        # 检查是否已处理
        output_file <- file.path(output_dir, paste0(sample_id, "_10000bp.segments.rds"))
        if (file.exists(output_file)) {
            return(list(status="already_processed", sample_id=sample_id))
        }

        # 创建进度文件
        progress_file <- file.path(progress_dir, paste0(sample_id, ".running"))
        file.create(progress_file)

        # 构建文件路径
        file_path <- file.path(input_dir,
                              paste0(sample_id, "_10000bp.correctedReadcount.rds"))

        # 检查文件是否存在
        if (!file.exists(file_path)) {
            unlink(progress_file)
            return(list(status="file_not_found", sample_id=sample_id))
        }

        # 读取数据
        corrected_reads <- readRDS(file_path)

        # 初始化结果列表
        all_segs <- list()

        # 处理每个染色体
        for(chr in c(1:22, "X")) {
            chr_name <- paste0("chr", chr)
            chr_data <- corrected_reads[corrected_reads$chr == chr_name, ]

            if(nrow(chr_data) > 0) {
                chr_data$chr <- gsub("chr", "", chr_data$chr)
                chr_data$chr <- as.factor(chr_data$chr)

                default_param <- HMMsegment(chr_data, getparam = TRUE)
                longseg_param <- default_param
                longseg_param$e <- 0.999999999999999
                longseg_param$strength <- 1e30
                longseg_segments <- HMMsegment(chr_data, longseg_param, verbose = FALSE)

                all_segs[[chr_name]] <- longseg_segments$segs
            }
        }

        # 保存结果
        saveRDS(all_segs, output_file)

        # 清理内存
        rm(corrected_reads, all_segs)
        gc()

        # 删除进度文件
        unlink(progress_file)

        print(sprintf("[%s] 成功处理样本: %s",
                     format(Sys.time(), "%H:%M:%S"),
                     sample_id))

        return(list(status="success", sample_id=sample_id))

    }, error = function(e) {
        print(sprintf("[%s] 处理样本 %s 时出错: %s",
                     format(Sys.time(), "%H:%M:%S"),
                     sample_id,
                     e$message))
        unlink(file.path(progress_dir, paste0(sample_id, ".running")))
        return(list(status="error", sample_id=sample_id, error=e$message))
    })
}

# 设置并行核心数
num_cores <- 4  # 固定使用4核
print(sprintf("使用核心数: %d", num_cores))

# 设置批处理
batch_size <- 8  # 每批8个样本
sample_batches <- split(remaining_samples,
                       ceiling(seq_along(remaining_samples)/batch_size))
total_batches <- length(sample_batches)
start_time <- Sys.time()

print(sprintf("\n总批次数: %d, 每批样本数: %d", total_batches, batch_size))

for(batch_idx in seq_along(sample_batches)) {
    batch <- sample_batches[[batch_idx]]
    current_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(current_time, start_time, units="secs"))

    # 计算预计剩余时间
    if(batch_idx > 1) {
        eta <- calculate_eta(batch_idx - 1, total_batches, elapsed_time)
    } else {
        eta <- "计算中..."
    }

    print(sprintf("\n========== 批次进度 %d/%d (%.1f%%) ==========",
                 batch_idx,
                 total_batches,
                 batch_idx/total_batches*100))
    print(sprintf("当前内存使用: %.2f GB", memory.size()/1024))
    print(sprintf("预计剩余时间: %s", eta))
    print(sprintf("当前批次样本数: %d", length(batch)))

    cl <- makeCluster(num_cores)
    clusterExport(cl, c("process_sample", "input_dir", "output_dir", "progress_dir"))
    clusterEvalQ(cl, {
        library(HMMcopy)
    })

    results_batch <- parLapply(cl, batch, process_sample)
    stopCluster(cl)

    # 输出当前批次处理结果统计
    success_count <- sum(sapply(results_batch, function(x) x$status == "success"))
    print(sprintf("\n当前批次处理完成:"))
    print(sprintf("成功: %d/%d", success_count, length(batch)))
    print(sprintf("失败: %d/%d", length(batch) - success_count, length(batch)))

    # 清理内存
    rm(results_batch)
    gc()
    Sys.sleep(5)
}

print(sprintf("\n========== 所有处理完成 =========="))
print(sprintf("总耗时: %.1f 分钟",
             as.numeric(difftime(Sys.time(), start_time, units="mins"))))

# 最终处理结果统计
final_processed_samples <- get_processed_samples()
final_remaining <- setdiff(all_samples, final_processed_samples)

print(sprintf("\n最终统计:"))
print(sprintf("总样本数: %d", length(all_samples)))
print(sprintf("成功处理: %d", length(final_processed_samples)))
print(sprintf("未处理/失败: %d", length(final_remaining)))

# 如果有未处理完的样本，保存到文件
if(length(final_remaining) > 0) {
    failed_samples <- data.frame(sample_id = final_remaining)
    write.csv(failed_samples,
              file.path(output_dir, "failed_samples.csv"),
              row.names=FALSE)
    print(sprintf("\n未处理完的样本已保存到: %s",
                 file.path(output_dir, "failed_samples.csv")))
}