
# 检查并安装必要的包
if (!requireNamespace("HMMcopy", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("HMMcopy")
}
library(HMMcopy)

# 定义窗口大小和对应的文件夹名称
window_sizes <- c(10000)
window_folders <- c("window_10000")

# 定义基础路径
base_path <- "/share/home/lsy_student/chenyanchao/GDM/CNV/HMMcopy/preload"

# 定义输出的基础路径
output_base_path <- "/share/home/lsy_chenyanchao/projects/hmmcopy/samples120k/hmmcopy_out/Longgang/hmmcopy_1"

# 定义各类文件的路径
gc_path <- file.path(base_path, "gc")
mappability_path <- file.path(base_path, "mappability")
readcounts_path <- "/share/home/lsy_chenyanchao/projects/hmmcopy/samples120k/readcounts_wigs/Longgang/Longgang_wigs_1"

# 定义函数来获取样本文件
get_sample_files <- function(window_folder) {
  pattern <- "\\.readcounts\\.wig$"
  files <- list.files(file.path(readcounts_path, window_folder), pattern = pattern, full.names = TRUE)

  # 输出找到的文件数目和示例
  cat("Found", length(files), "sample files in", window_folder, "\n")
  if (length(files) > 0) {
    cat("Sample file example:", files[1], "\n")
  }

  return(files)
}

# 主处理函数
process_sample <- function(sample_file, window_size, window_folder) {
  # 构建输出文件路径（移到函数开头以便提前检查）
  filename <- basename(sample_file)  # 提取文件名
  sample <- sub("\\.readcounts\\.wig$", "", filename)  # 去除后缀
  output_dir_sample <- file.path(output_base_path, window_folder)
  output_file <- file.path(output_dir_sample, paste0(sample, ".correctedReadcount.rds"))

  # 检查输出文件是否已存在
  if (file.exists(output_file)) {
    cat("跳过已处理的样本:", sample, "\n")
    return(NULL)
  }

  # 构建gc和mappability文件路径
  gc_file <- file.path(gc_path, window_folder, paste0("Homo_sapiens_assembly38_", window_size, "bp.gc_filtered.wig"))
  map_file <- file.path(mappability_path, window_folder, paste0("k100.Umap.MultiTrackMappability_", window_size, "bp.wig"))

  # 输出构建的文件路径以便调试
  cat("Constructed GC file path:", gc_file, "\n")
  cat("Constructed Mappability file path:", map_file, "\n")

  # 检查文件是否存在
  if (!file.exists(sample_file)) {
    cat("Error: Sample file not found:", sample_file, "\n")
    return(NULL)
  }
  if (!file.exists(gc_file)) {
    cat("Error: GC file not found:", gc_file, "\n")
    return(NULL)
  }
  if (!file.exists(map_file)) {
    cat("Error: Mappability file not found:", map_file, "\n")
    return(NULL)
  }

  # 读取WIG文件
  cat("Reading readcounts from:", sample_file, "\n")
  readcount <- tryCatch({
    HMMcopy::wigsToRangedData(sample_file, gc_file, map_file)
  }, error = function(e) {
    cat("Error reading WIG files for", sample_file, ":", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(readcount)) {
    return(NULL)
  }

  # 从文件名中提取样本ID
  cat("Extracted sample ID:", sample, "\n")

  # 创建HMMcopy对象
  cat("Correcting readcount for sample:", sample, "\n")
  copy <- tryCatch({
    HMMcopy::correctReadcount(readcount)
  }, error = function(e) {
    cat("Error correcting readcounts for", sample, ":", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(copy)) {
    return(NULL)
  }

  # 定义输出文件路径
  if (!dir.exists(output_dir_sample)) {
    dir.create(output_dir_sample, recursive = TRUE)
    cat("Created output directory:", output_dir_sample, "\n")
  }

  # 保存 corrected readcount
  cat("Saving corrected readcount to:", output_file, "\n")
  tryCatch({
    saveRDS(copy, file = output_file)
  }, error = function(e) {
    cat("Error saving corrected readcount for", sample, ":", conditionMessage(e), "\n")
  })

  return(copy)
}

# 主循环
for (i in seq_along(window_sizes)) {
  window_size <- window_sizes[i]
  window_folder <- window_folders[i]

  cat("\n==============================\n")
  cat("处理窗口大小:", window_size, "来自文件夹:", window_folder, "\n")
  cat("==============================\n")

  sample_files <- get_sample_files(window_folder)

  if (length(sample_files) == 0) {
    cat("警告: 在", window_folder, "中未找到样本文件\n")
    next
  }

  # 添加进度计数
  total_samples <- length(sample_files)
  processed_samples <- 0

  for (sample_file in sample_files) {
    processed_samples <- processed_samples + 1
    cat("\n--------------------------------\n")
    cat(sprintf("处理进度: %d/%d (%d%%)\n",
                processed_samples, total_samples,
                round(processed_samples/total_samples*100)))
    cat("正在处理:", sample_file, "窗口大小:", window_size, "\n")
    cat("--------------------------------\n")

    tryCatch({
      segs <- process_sample(sample_file, window_size, window_folder)
    }, error = function(e) {
      cat("处理错误", sample_file, ":", conditionMessage(e), "\n")
    })
  }
}

cat("\n所有样本处理完成。\n")





