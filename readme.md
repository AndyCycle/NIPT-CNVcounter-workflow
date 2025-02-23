# CNV 矫正估计项目

本项目提供了一系列利用 **HMMcopy** 包对高通量测序数据进行拷贝数变异（CNV）矫正估计的脚本。这些脚本构成了一个完整的工作流，从 reads 计数、数据校正、拷贝数估计，到最终的矩阵生成，均以 R 脚本和 Shell 脚本的形式实现。项目的各个模块能够灵活配置参数，方便用户根据实际需求定制数据处理流程。

---

## 项目模块概览

- **readcounter_parallel.sh**
  利用 HMMcopy_utils 工具中的 `readCounter` 对 BAM 文件进行 reads 计数，生成无重叠窗口区间下的 `.wig` 文件。该脚本采用 GNU Parallel 以并行方式提交任务，并利用 SLURM 进行作业调度。
  **自定义设置部分：**
  - `your_conda_path`：Conda 环境路径。
  - `your_Read_Counter_path`：readCounter 可执行程序的完整路径。
  - `your_bam_dir`：存放 BAM 文件的根目录。
  - `your_output_dir`：输出 wig 文件的目录。
  - `your_window_size`：指定的窗口大小。

- **hmmcopy_readcounts.r**
  使用 **HMMcopy** 包读取已生成的 `.wig` 文件，对 readcounts 数据进行校正，输出 `.rds` 格式的校正结果。
  **自定义设置部分：**
  - `base_path`：原始数据文件的根目录。
  - `output_base_path`：校正后的数据保存目录。
  - `readcounts_path`：存放 `.wig` 文件的文件夹路径。
  - `window_sizes` 与 `window_folders`：窗口大小及其对应的文件夹名称。

- **hmmcopy_copy.r**
  基于校正后的 readcounts 数据，通过调用 **HMMsegment** 函数，进一步计算各区间的拷贝数（输出中 state 列数据，实际拷贝数为 state - 1），结果保存为 `.rds` 格式。
  **自定义设置部分：**
  - 设置输入和输出目录（对应校正数据与拷贝数结果）。
  - 修改样本文件名匹配规则和进程管理（并行处理和任务进度控制）。

- **readcounts_table.r**
  整合所有校正后的 readcounts 数据，生成包含各样本 readcounts 信息的矩阵。该矩阵以行保存窗口（包含染色体、起始和终止位置），每个样本的数据作为一列存入。
  **自定义设置部分：**
  - `base_input_dir` 与 `base_output_dir`：输入、输出基础目录。
  - `hospital_name`、`your_window_size`、`your_hmmcopy_number`：医院名称、窗口大小和对应的 hmmcopy 序号。
  - 代码内部使用 `file.path()` 构建了标准目录路径，便于跨平台使用。

- **copy_table.r**
  利用预处理得到的包含拷贝数信息的 `.segments.rds` 文件，生成最终的 CNV 矩阵（cnv_table.rds）。
  **自定义设置部分：**
  - 输入参数包括 `base_input_dir`、`hospital_name`、`window_size`、`segs_number` 等。
  - 构建输入输出路径时，根据参数生成对应的分段和窗口目录。
  - 支持断点续处理，通过记录已处理样本文件来实现任务进度维护。

- **generate_scripts.r**
  一个生成器脚本，用于自动生成多个任务脚本与相应的 sbatch 作业提交脚本。
  用户可通过设定医院名称、窗口大小和分段数等参数，快速生成分布在不同文件夹下的任务脚本，方便在集群上批量提交。
  **自定义设置部分：**
  - 修改 `hospital_names`、`window_sizes`、`segs_count` 参数。
  - 设置工作目录（`setwd("your_working_directory")`）以及生成脚本存放的目录。
  - 模板中的文件路径（如 `your_base_input_dir`、`your_output_base_dir`、`your_conda_path`、`your_working_directory` 等）需根据实际环境进行调整。

---

## 文件目录结构示例

项目文件的存储方案如下，供参考和部署时调整：
