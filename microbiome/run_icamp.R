# 进入 iCAMP 目录
setwd("/bmp/backup/zhaosy/ws/china_16s_pipeline/scripts/iCAMP/R")

# 加载参数
source("metaset.R")

# 加载 iCAMP 主脚本
# source("icamp.cm.r")
source("icamp.big.r")

# 读入数据
comm <- read.table(
  "/bmp/backup/zhaosy/ws/china_16s_pipeline/results/diversity_analysis/betaNTI/pa_otu_table.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

tree <- ape::read.tree("/bmp/backup/zhaosy/ws/china_16s_pipeline/results/feature_table_from_Jiahe/asv_fasttree.nwk")

# 运行 iCAMP
# res <- icamp.cm(
#   comm = comm,
#   tree = tree,
#   prefix = "iCAMP_result",
#   bin.size = 24,     # 推荐 20–24
#   rand.time = rand.time
# )
res <- icamp.big(
  comm = comm,
  tree = tree,
  prefix = "iCAMP_result",
  bin.size = 24,     # 推荐 20–24
)

# 查看最终五过程结果
print(res$ProcessImportance)


# 设置输出目录
outdir <- "/bmp/backup/zhaosy/ws/china_16s_pipeline/results/diversity_analysis/betaNTI/iCAMP_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 保存五过程
write.csv(
  res$ProcessImportance,
  file = file.path(outdir, "ProcessImportance.csv"),
  quote = FALSE,
  row.names = TRUE
)

# 保存 betaNTI 样本对结果
write.csv(
  res$BN.pairwise.result,
  file = file.path(outdir, "betaNTI_pairwise.csv"),
  quote = FALSE,
  row.names = FALSE
)

# 保存 betaNTI bin-level结果
write.csv(
  res$BN.bin.result,
  file = file.path(outdir, "betaNTI_bin.csv"),
  quote = FALSE,
  row.names = FALSE
)

# 保存 RCbray bin-level结果
if(!is.null(res$RC.bin.result)){
  write.csv(
    res$RC.bin.result,
    file = file.path(outdir, "RC_bin.csv"),
    quote = FALSE,
    row.names = FALSE
  )
}

# 保存整个 res 对象方便后续分析
saveRDS(res, file = file.path(outdir, "iCAMP_res.rds"))

# 打印完成信息
cat("所有 iCAMP 结果已保存到：", outdir, "\n")
