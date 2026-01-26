required_pkgs <- c("picante", "ape")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(picante)
library(ape)

calculate_betaNTI <- function(
  otu_file,
  tree_file,
  nreps = 999,
  abundance_weighted = FALSE,
  null_model = "taxa.labels",
  output_file = NULL
) {

  # --------------------------------------------------
  # 1. 读取数据
  # --------------------------------------------------
  comm <- read.table(
    otu_file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )

  tree <- read.tree(tree_file)

  # --------------------------------------------------
  # 2. 数据校验
  # --------------------------------------------------
  taxa_comm <- colnames(comm)
  taxa_tree <- tree$tip.label

  if (!setequal(taxa_comm, taxa_tree)) {
    stop(
      "Taxa mismatch between OTU table and tree\n",
      "Only in OTU: ", paste(setdiff(taxa_comm, taxa_tree), collapse = ", "), "\n",
      "Only in tree: ", paste(setdiff(taxa_tree, taxa_comm), collapse = ", ")
    )
  }

  # 按 tree 顺序重排 OTU 表
  comm <- comm[, tree$tip.label]

  # --------------------------------------------------
  # 3. betaMNTD / betaNTI 计算
  # --------------------------------------------------
  phylo_dist <- cophenetic(tree)

  result <- ses.beta.mntd(
    samp = comm,
    dis = phylo_dist,
    null.model = null_model,
    abundance.weighted = abundance_weighted,
    runs = nreps
  )

  # --------------------------------------------------
  # 4. 生态过程分类（Stegen）
  # --------------------------------------------------
  result$process <- ifelse(
    result$mntd.obs.z > 2, "Variable selection",
    ifelse(
      result$mntd.obs.z < -2, "Homogeneous selection",
      "Stochastic"
    )
  )

  # --------------------------------------------------
  # 5. 输出
  # --------------------------------------------------
  if (!is.null(output_file)) {
    write.csv(result, output_file, quote = FALSE)
  }

  return(result)
}


if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    stop("Usage: Rscript betaNTI_module.R otu_table.tsv tree.nwk [output.csv]")
  }

  otu_file <- args[1]
  tree_file <- args[2]
  output_file <- ifelse(length(args) >= 3, args[3], "betaNTI_result.csv")

  calculate_betaNTI(
    otu_file = otu_file,
    tree_file = tree_file,
    output_file = output_file
  )
}

