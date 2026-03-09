#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
})



# ========== 参数 ==========
option_list <- list(
  make_option("--comm", type="character", help="community matrix csv"),
  make_option("--env", type="character", help="environment data csv"),
  make_option("--formula", type="character", help="env variables, e.g. pH + Temp"),
  make_option("--distance", type="character", default="bray"),
  make_option("--outdir", type="character", default="dbRDA_results")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ========== 读取数据 ==========
# comm <- read.csv(opt$comm, row.names = 1, check.names = FALSE)
# env  <- read.csv(opt$env,  row.names = 1, check.names = FALSE)
comm <- read.table(
  opt$comm,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

env <- read.table(
  opt$env,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# 保证样本一致
comm <- comm[rownames(env), ]

# 去除全 0 物种
comm <- comm[, colSums(comm) > 0]

# Hellinger 转换（推荐）
comm_hel <- decostand(comm, method = "hellinger")

# ========== db-RDA ==========
form <- as.formula(paste("comm_hel ~", opt$formula))

db <- capscale(
  form,
  data = env,
  distance = opt$distance
)

# ========== 统计检验 ==========
anova_all   <- anova(db, permutations = 999)
anova_terms <- anova(db, by = "terms", permutations = 999)
anova_axis  <- anova(db, by = "axis", permutations = 999)

# ========== 提取结果 ==========
# 样本坐标
sites <- scores(db, display = "sites", scaling = 2)
sites <- as.data.frame(sites)
sites$SampleID <- rownames(sites)

# 环境变量箭头
bp <- scores(db, display = "bp", scaling = 2)
bp <- as.data.frame(bp)
bp$Variable <- rownames(bp)

# 轴解释率
eig <- summary(db)$cont$importance
eig <- as.data.frame(t(eig))

# ========== 保存 ==========
write.csv(sites, file.path(opt$outdir, "sites_scores.csv"), row.names = FALSE)
write.csv(bp,    file.path(opt$outdir, "env_vectors.csv"), row.names = FALSE)
write.csv(eig,   file.path(opt$outdir, "axis_importance.csv"))

write.csv(as.data.frame(anova_all),
          file.path(opt$outdir, "anova_overall.csv"))
write.csv(as.data.frame(anova_terms),
          file.path(opt$outdir, "anova_terms.csv"))
write.csv(as.data.frame(anova_axis),
          file.path(opt$outdir, "anova_axis.csv"))

cat("db-RDA analysis finished.\n")
