library(HMMcopy)
library(data.table)
library(parallel)

# 设置输入输出目录
input_dir <- "/share/home/lsy_chenyanchao/projects/hmmcopy/samples120k/hmmcopy_out/Baoan/hmmcopy_1/window_10000"
output_dir <- "/share/home/lsy_chenyanchao/projects/hmmcopy/samples120k/hmmcopy_segs/Baoan/segs_1/window_10000"
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

# 获取待处理样本
all_samples <- get_all_samples()
processed_samples <- get_processed_samples()
remaining_samples <- setdiff(all_samples, processed_samples)

print(sprintf("Total input samples: %d", length(all_samples)))
print(sprintf("Already processed samples: %d", length(processed_samples)))
print(sprintf("Remaining samples to process: %d", length(remaining_samples)))

# 定义处理单个样本的函数
process_sample <- function(sample_id) {
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

        # 删除进度文件
        unlink(progress_file)

        return(list(status="success", sample_id=sample_id))

    }, error = function(e) {
        unlink(file.path("progress", paste0(sample_id, ".running")))
        return(list(status="error", sample_id=sample_id, error=e$message))
    })
}

# 设置并行核心数（根据服务器配置调整）
num_cores <- min(8, detectCores() - 2)  # 留出2个核心给系统
print(sprintf("Using %d cores for parallel processing", num_cores))

# 将样本分批处理
batch_size <- 16  # 每批处理16个样本
sample_batches <- split(remaining_samples,
                       ceiling(seq_along(remaining_samples)/batch_size))

for(batch in sample_batches) {
    print(sprintf("Processing batch with %d samples", length(batch)))

    cl <- makeCluster(num_cores)
    clusterExport(cl, c("process_sample", "input_dir", "output_dir", "progress_dir"))
    clusterEvalQ(cl, {
        library(HMMcopy)
    })

    results_batch <- parLapply(cl, batch, process_sample)
    stopCluster(cl)

    # 处理完一批后强制清理内存
    gc()
}

# 执行并行处理
results <- parLapply(cl, remaining_samples, process_sample)
stopCluster(cl)

# 处理结果统计
process_results <- data.frame(
    sample_id = sapply(results, function(x) x$sample_id),
    status = sapply(results, function(x) x$status)
)

# 保存处理结果统计
write.csv(process_results, "processing_summary.csv", row.names=FALSE)

# 检查是否有处理失败的样本
failed_samples <- process_results[process_results$status != "success", ]
if(nrow(failed_samples) > 0) {
    print("Some samples failed to process:")
    print(failed_samples)
    write.csv(failed_samples, file.path(output_dir, "failed_samples.csv"), row.names=FALSE)
}

# 最终统计
print(sprintf("Processing completed:"))
print(sprintf("Total samples: %d", nrow(process_results)))
print(sprintf("Successful: %d", sum(process_results$status == "success")))
print(sprintf("Failed: %d", sum(process_results$status != "success")))

# 保存处理结果统计
write.csv(process_results, file.path(output_dir, "processing_summary.csv"), row.names=FALSE)
