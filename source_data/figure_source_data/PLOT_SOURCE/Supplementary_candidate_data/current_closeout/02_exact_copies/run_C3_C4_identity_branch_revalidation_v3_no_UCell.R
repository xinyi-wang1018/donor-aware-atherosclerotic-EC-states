#!/usr/bin/env Rscript

# ============================================================================
# 动脉粥样硬化 EC：C3/C4 血管身份与 C7→C3→双分支重新验证
# 输入对象：/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds
# 输出目录：/root/AS_EC_identity_branch_revalidation/
#
# 说明：身份程序评分使用 Seurat::AddModuleScore，不再依赖 UCell。
#
# 主要完成：
# 1) C3 是否仍保留动脉来源，还是毛细血管/血管生成型 EC；
# 2) C4 是否是保留内皮身份的 ECM/partial EndMT-like，还是静脉型 EC/污染；
# 3) 在未校正 RNA-PCA 与 Harmony 空间分别重建局部结构；
# 4) 用 Slingshot 在不指定终点的情况下检查两条分支是否自然恢复；
# 5) leave-one-dataset-out 与患者/样本级 bootstrap 检查分支稳定性。
# ============================================================================

options(stringsAsFactors = FALSE, timeout = 3600, download.file.method = "libcurl")
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
set.seed(20260728)

# ------------------------------ 0. 参数 -------------------------------------
INPUT_RDS <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
OUTDIR   <- "/root/AS_EC_identity_branch_revalidation"

TARGET_CLUSTERS <- c("1", "2", "3", "4", "5", "7")
DISPLAY_ORDER   <- c("7", "3", "5", "1", "2", "4")
ROOT_CLUSTER    <- "7"

N_HVG        <- 3000
N_PCS        <- 40
DIMS_USE     <- 1:30
RESOLUTION   <- 0.60
N_BOOT       <- 50       # 患者/样本级 bootstrap 次数
BOOT_FRAC    <- 0.80     # 每次抽取 80% 样本
RUN_LODO     <- TRUE     # 5 次 leave-one-dataset-out
RUN_BOOTSTRAP <- TRUE
INSTALL_MISSING <- FALSE
N_CORES <- max(1L, min(8L, parallel::detectCores(logical = FALSE)))

# -------------------------- 1. 创建目录与日志 -------------------------------
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
subdirs <- c(
  "01_identity_scores", "02_marker_coexpression", "03_composition",
  "04_local_reanalysis", "05_slingshot", "06_LODO", "07_bootstrap",
  "08_final_summary"
)
invisible(lapply(file.path(OUTDIR, subdirs), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

LOG_FILE <- file.path(OUTDIR, "00_run_log.txt")
zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

on.exit({
  cat("\n===== sessionInfo =====\n")
  print(sessionInfo())
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("Start time:", format(Sys.time()), "\n")
cat("Input:", INPUT_RDS, "\n")
cat("Output:", OUTDIR, "\n")
cat("Cores:", N_CORES, "\n\n")

# ---------------------------- 2. 安装/加载包 --------------------------------
cran_pkgs <- c(
  "Seurat", "ggplot2", "dplyr", "tidyr", "patchwork", "pheatmap",
  "data.table", "harmony", "RColorBrewer", "scales"
)
bioc_pkgs <- c("slingshot", "SingleCellExperiment")

install_cran_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    if (!INSTALL_MISSING) stop("Missing CRAN packages: ", paste(miss, collapse = ", "))
    install.packages(miss, repos = "https://cloud.r-project.org", Ncpus = N_CORES)
  }
}

install_bioc_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    if (!INSTALL_MISSING) stop("Missing Bioconductor packages: ", paste(miss, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(miss, ask = FALSE, update = FALSE, Ncpus = N_CORES)
  }
}

install_cran_if_missing(cran_pkgs)
install_bioc_if_missing(bioc_pkgs)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
  library(data.table)
  library(harmony)
  library(slingshot)
  library(SingleCellExperiment)
  library(Matrix)
})

# ----------------------------- 3. 通用函数 -----------------------------------
first_present <- function(candidates, available, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% available]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "))
  }
  return(NA_character_)
}

normalize_cluster_label <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub("^cluster[_ -]*", "", y, ignore.case = TRUE)
  y <- gsub("^[Cc]", "", y)
  num <- regmatches(y, regexpr("[0-9]+$", y))
  use_num <- nchar(num) > 0
  y[use_num] <- num[use_num]
  y
}

join_rna_layers_if_needed <- function(obj) {
  if (!"RNA" %in% names(obj@assays)) stop("RNA assay not found.")
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(
    SeuratObject::JoinLayers(obj, assay = "RNA"),
    error = function(e) {
      message("JoinLayers skipped: ", conditionMessage(e))
      obj
    }
  )
  obj
}

get_assay_matrix <- function(obj, assay = "RNA", layer = c("data", "counts")) {
  layer <- match.arg(layer)
  tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e1) {
      tryCatch(
        SeuratObject::GetAssayData(obj, assay = assay, slot = layer),
        error = function(e2) stop("Cannot obtain ", assay, " ", layer, " matrix: ", conditionMessage(e2))
      )
    }
  )
}

safe_ggsave <- function(filename, plot, width, height, dpi = 300) {
  ggsave(filename = filename, plot = plot, width = width, height = height,
         dpi = dpi, bg = "white", limitsize = FALSE)
}

save_csv <- function(x, path) {
  data.table::fwrite(as.data.frame(x), path)
}

contains_in_order <- function(sequence, required) {
  pos <- match(required, sequence)
  all(!is.na(pos)) && all(diff(pos) > 0)
}

lineage_summary <- function(sds) {
  lin <- slingshot::slingLineages(sds)
  if (length(lin) == 0) {
    return(data.frame(
      lineage = character(), sequence = character(), endpoint = character(),
      has_L1_exact = logical(), has_L1_core = logical(), has_L2 = logical()
    ))
  }
  out <- lapply(seq_along(lin), function(i) {
    seq_i <- as.character(lin[[i]])
    data.frame(
      lineage = paste0("Lineage", i),
      sequence = paste(seq_i, collapse = "->"),
      endpoint = tail(seq_i, 1),
      has_L1_exact = contains_in_order(seq_i, c("7", "3", "5", "1", "2")),
      has_L1_core  = contains_in_order(seq_i, c("7", "3", "2")),
      has_L2       = contains_in_order(seq_i, c("7", "3", "4")),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

run_slingshot_matrix <- function(embedding, labels, start_cluster = "7") {
  keep <- complete.cases(embedding) & !is.na(labels)
  emb <- as.matrix(embedding[keep, , drop = FALSE])
  lab <- factor(as.character(labels[keep]))
  if (!start_cluster %in% levels(lab)) {
    stop("Root cluster ", start_cluster, " not present.")
  }
  if (length(levels(droplevels(lab))) < 3) {
    stop("Fewer than 3 clusters remain; cannot infer branching.")
  }
  slingshot::slingshot(
    emb,
    clusterLabels = lab,
    start.clus = start_cluster,
    stretch = 0,
    reweight = TRUE,
    reassign = TRUE
  )
}

plot_slingshot_2d <- function(embedding, labels, sds, path, title_text) {
  labs <- factor(as.character(labels))
  cols <- setNames(grDevices::hcl.colors(length(levels(labs)), "Dark 3"), levels(labs))
  png(path, width = 2400, height = 2000, res = 250)
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(
    embedding[, 1], embedding[, 2],
    col = adjustcolor(cols[as.character(labs)], alpha.f = 0.55),
    pch = 16, cex = 0.55,
    xlab = colnames(embedding)[1], ylab = colnames(embedding)[2],
    main = title_text
  )
  try(lines(sds, lwd = 4, col = "black"), silent = TRUE)
  legend("topright", legend = names(cols), col = cols, pch = 16,
         title = "EC cluster", cex = 0.85, bty = "n")
  dev.off()
}

add_pseudotime_and_plot <- function(obj, sds, umap_reduction, prefix, outdir) {
  pt <- slingshot::slingPseudotime(sds, na = TRUE)
  if (is.null(dim(pt))) pt <- matrix(pt, ncol = 1)
  rownames(pt) <- rownames(Embeddings(obj, umap_reduction))
  for (j in seq_len(ncol(pt))) {
    nm <- paste0(prefix, "_PT", j)
    obj[[nm]] <- pt[Cells(obj), j]
    um <- Embeddings(obj, umap_reduction)
    df <- data.frame(
      x = um[, 1], y = um[, 2], pseudotime = pt[rownames(um), j]
    )
    p <- ggplot(df, aes(x, y, color = pseudotime)) +
      geom_point(size = 0.55, alpha = 0.8) +
      scale_color_viridis_c(na.value = "grey85") +
      theme_classic(base_size = 13) +
      labs(title = paste0(prefix, " lineage ", j), x = "UMAP_1", y = "UMAP_2")
    safe_ggsave(file.path(outdir, paste0(prefix, "_lineage", j, "_pseudotime_UMAP.png")),
                p, 7.5, 6.5)
  }
  obj
}

run_harmony_robust <- function(obj, dataset_col, dims_use, reduction_save = "harmony_local") {
  if (is.na(dataset_col) || length(unique(obj[[dataset_col]][, 1])) < 2) {
    message("Harmony skipped: fewer than two dataset levels.")
    return(obj)
  }
  obj <- tryCatch(
    harmony::RunHarmony(
      object = obj,
      group.by.vars = dataset_col,
      reduction.use = "pca",
      dims.use = dims_use,
      reduction.save = reduction_save,
      project.dim = FALSE,
      verbose = TRUE
    ),
    error = function(e_old) {
      message("Old Harmony interface failed; retrying new interface: ", conditionMessage(e_old))
      harmony::RunHarmony(
        object = obj,
        group.by.vars = dataset_col,
        reduction = "pca",
        dims.use = dims_use,
        reduction.save = reduction_save,
        project.dim = FALSE,
        verbose = TRUE
      )
    }
  )
  obj
}

prepare_local_object <- function(obj, dataset_col, run_harmony = TRUE) {
  DefaultAssay(obj) <- "RNA"
  obj <- join_rna_layers_if_needed(obj)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                       verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_HVG,
                              verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  max_pcs <- min(N_PCS, length(VariableFeatures(obj)) - 1L, ncol(obj) - 1L)
  if (max_pcs < 10) stop("Too few cells/features for PCA.")
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = max_pcs,
                verbose = FALSE, seed.use = 20260728)
  dims_use <- seq_len(min(max(DIMS_USE), max_pcs))

  obj <- RunUMAP(
    obj, reduction = "pca", dims = dims_use,
    reduction.name = "umap_pca_local", reduction.key = "UPCA_",
    n.neighbors = 30, min.dist = 0.30, metric = "cosine",
    seed.use = 20260728, verbose = FALSE
  )
  obj <- FindNeighbors(
    obj, reduction = "pca", dims = dims_use,
    graph.name = c("pca_nn", "pca_snn"), verbose = FALSE
  )
  old_ident <- Idents(obj)
  obj <- FindClusters(obj, graph.name = "pca_snn", resolution = RESOLUTION,
                      algorithm = 1, random.seed = 20260728, verbose = FALSE)
  obj$local_cluster_pca <- as.character(Idents(obj))
  Idents(obj) <- old_ident

  if (run_harmony && !is.na(dataset_col) && length(unique(obj[[dataset_col]][, 1])) >= 2) {
    obj <- run_harmony_robust(obj, dataset_col, dims_use, reduction_save = "harmony_local")
    obj <- RunUMAP(
      obj, reduction = "harmony_local", dims = dims_use,
      reduction.name = "umap_harmony_local", reduction.key = "UHAR_",
      n.neighbors = 30, min.dist = 0.30, metric = "cosine",
      seed.use = 20260728, verbose = FALSE
    )
    obj <- FindNeighbors(
      obj, reduction = "harmony_local", dims = dims_use,
      graph.name = c("harmony_nn", "harmony_snn"), verbose = FALSE
    )
    old_ident <- Idents(obj)
    obj <- FindClusters(obj, graph.name = "harmony_snn", resolution = RESOLUTION,
                        algorithm = 1, random.seed = 20260728, verbose = FALSE)
    obj$local_cluster_harmony <- as.character(Idents(obj))
    Idents(obj) <- old_ident
  }
  obj
}

aggregate_graph_by_cluster <- function(obj, graph_name, cluster_col) {
  if (!graph_name %in% names(obj@graphs)) return(NULL)
  g <- obj@graphs[[graph_name]]
  cells <- intersect(rownames(g), Cells(obj))
  g <- g[cells, cells, drop = FALSE]
  labs <- factor(as.character(obj[[cluster_col]][cells, 1]))
  mm <- Matrix::sparse.model.matrix(~ 0 + labs)
  colnames(mm) <- sub("^labs", "", colnames(mm))
  agg <- t(mm) %*% g %*% mm
  agg <- as.matrix(agg)
  row_prop <- agg / pmax(rowSums(agg), 1e-12)
  list(raw = agg, row_prop = row_prop)
}

trajectory_flags <- function(summary_df) {
  data.frame(
    n_lineages = nrow(summary_df),
    L1_exact = any(summary_df$has_L1_exact),
    L1_core = any(summary_df$has_L1_core),
    L2 = any(summary_df$has_L2),
    both_L1core_L2 = any(summary_df$has_L1_core) && any(summary_df$has_L2),
    stringsAsFactors = FALSE
  )
}

# ------------------------------ 4. 读取对象 ----------------------------------
if (!file.exists(INPUT_RDS)) stop("Input RDS not found: ", INPUT_RDS)
obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.")
obj <- join_rna_layers_if_needed(obj)
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = FALSE)

meta_cols <- colnames(obj@meta.data)
cluster_col <- first_present(
  c("EC_cluster", "ec_cluster", "EC_clusters", "seurat_clusters"),
  meta_cols, required = TRUE, label = "EC cluster column"
)
dataset_col <- first_present(
  c("dataset_id", "dataset", "study", "GSE", "batch", "source_dataset"),
  meta_cols, required = FALSE, label = "dataset column"
)
sample_col <- first_present(
  c("sample_id", "sample", "orig.ident", "patient_id", "donor_id", "subject_id"),
  meta_cols, required = FALSE, label = "sample column"
)
tissue_col <- first_present(
  c("tissue", "tissue_type", "vascular_bed", "artery_site", "site"),
  meta_cols, required = FALSE, label = "tissue column"
)
disease_col <- first_present(
  c("disease", "condition", "group", "disease_status", "phenotype", "status"),
  meta_cols, required = FALSE, label = "disease column"
)

obj$EC_cluster_revalidated <- normalize_cluster_label(obj[[cluster_col]][, 1])
cluster_col_use <- "EC_cluster_revalidated"

present_clusters <- sort(unique(obj[[cluster_col_use]][, 1]))
missing_target <- setdiff(TARGET_CLUSTERS, present_clusters)
if (length(missing_target) > 0) {
  stop("Target clusters absent from object: ", paste(missing_target, collapse = ", "),
       "\nPresent clusters: ", paste(present_clusters, collapse = ", "))
}

metadata_report <- c(
  paste0("input_rds=", INPUT_RDS),
  paste0("n_cells=", ncol(obj)),
  paste0("n_genes=", nrow(obj)),
  paste0("original_cluster_col=", cluster_col),
  paste0("working_cluster_col=", cluster_col_use),
  paste0("dataset_col=", ifelse(is.na(dataset_col), "NOT_FOUND", dataset_col)),
  paste0("sample_col=", ifelse(is.na(sample_col), "NOT_FOUND", sample_col)),
  paste0("tissue_col=", ifelse(is.na(tissue_col), "NOT_FOUND", tissue_col)),
  paste0("disease_col=", ifelse(is.na(disease_col), "NOT_FOUND", disease_col)),
  paste0("present_clusters=", paste(present_clusters, collapse = ","))
)
writeLines(metadata_report, file.path(OUTDIR, "00_detected_metadata_columns.txt"))
cat(paste(metadata_report, collapse = "\n"), "\n\n")

# 只保留主轨迹候选 C1/C2/C3/C4/C5/C7
cells_target <- Cells(obj)[obj[[cluster_col_use]][, 1] %in% TARGET_CLUSTERS]
target <- subset(obj, cells = cells_target)
target[[cluster_col_use]] <- factor(target[[cluster_col_use]][, 1], levels = DISPLAY_ORDER)
Idents(target) <- cluster_col_use

# ------------------------- 5. 血管身份基因程序评分 ---------------------------
# 注意：静脉核心与静脉修复分开，避免仅凭 POSTN/MAFB 就循环论证 C4 是静脉型。
signatures <- list(
  Pan_endothelial = c("PECAM1", "VWF", "KDR", "CLDN5", "CDH5", "ESAM", "ENG", "RAMP2", "EGFL7", "ERG"),
  Arterial_core = c("SOX17", "EFNB2", "GJA4", "GJA5", "HEY1", "DLL4", "BMX", "CXCL12", "FBLN5", "SEMA3G"),
  Capillary_core = c("CA4", "RGCC", "EMCN", "PLVAP", "EDNRB", "BTNL9", "SLC9A3R2", "KDR", "FLT1"),
  Angiogenic_activation = c("APLN", "PGF", "ANGPT2", "ESM1", "ADM", "KDR", "FLT1", "CXCR4", "DLL4", "VEGFC"),
  Venous_core = c("NR2F2", "EPHB4", "ACKR1", "PLAT", "AQP1", "NRP2", "VWA1", "ADGRG6"),
  Venous_repair_like = c("NR2F2", "MAFB", "POSTN", "PLAT", "COL15A1", "VWA1", "AQP1"),
  ECM_remodeling = c("FN1", "BGN", "SPARC", "POSTN", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL8A1", "ELN"),
  EndMT_regulatory = c("VIM", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "TAGLN", "ACTA2", "COL8A1", "BGN"),
  Fibroblast = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "PDGFRA", "CFD", "C7", "COL6A1", "COL6A2"),
  SMC_pericyte = c("ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "CSPG4", "NOTCH3", "MCAM", "DES"),
  Inflammatory_activation = c("VCAM1", "SELE", "ICAM1", "CXCL2", "CXCL8", "CCL2", "NFKBIA", "STAT1", "IRF1")
)

all_genes <- rownames(target)
signatures_present <- lapply(signatures, function(gs) intersect(gs, all_genes))
sig_availability <- bind_rows(lapply(names(signatures), function(nm) {
  data.frame(
    signature = nm,
    requested_genes = length(signatures[[nm]]),
    present_genes = length(signatures_present[[nm]]),
    genes_present = paste(signatures_present[[nm]], collapse = ";"),
    genes_missing = paste(setdiff(signatures[[nm]], all_genes), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
save_csv(sig_availability,
         file.path(OUTDIR, "01_identity_scores", "01_signature_gene_availability.csv"))

valid_signatures <- signatures_present[vapply(signatures_present, length, integer(1)) >= 3]
if (length(valid_signatures) < 5) stop("Too few valid signatures after gene matching.")

# 不再依赖 UCell：直接使用 Seurat 自带的 AddModuleScore。
# 每个基因程序的分数 = 该程序基因平均表达 - 匹配表达水平的对照基因平均表达。
# 这样可以绕过 UCell 下载失败，同时仍可用于 C3/C4 身份程序的相对比较。
set.seed(20260728)
old_meta_cols <- colnames(target@meta.data)
target <- Seurat::AddModuleScore(
  object = target,
  features = unname(valid_signatures),
  assay = "RNA",
  name = "IdentityModule",
  nbin = 24,
  ctrl = 100,
  seed = 20260728,
  search = FALSE
)
raw_score_cols <- setdiff(colnames(target@meta.data), old_meta_cols)
expected_raw <- paste0("IdentityModule", seq_along(valid_signatures))
if (all(expected_raw %in% colnames(target@meta.data))) {
  raw_score_cols <- expected_raw
}
if (length(raw_score_cols) != length(valid_signatures)) {
  stop(
    "AddModuleScore returned an unexpected number of score columns. Expected ",
    length(valid_signatures), ", got ", length(raw_score_cols), ". Columns: ",
    paste(raw_score_cols, collapse = ", ")
  )
}
score_cols <- paste0(names(valid_signatures), "_ModuleScore")
for (i in seq_along(raw_score_cols)) {
  target[[score_cols[i]]] <- target[[raw_score_cols[i]]][, 1]
}
target@meta.data <- target@meta.data[, setdiff(colnames(target@meta.data), raw_score_cols), drop = FALSE]

# 5.1 小提琴图
score_df <- target@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, all_of(cluster_col_use), all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols), names_to = "signature", values_to = "score")

p_violin <- ggplot(score_df, aes(x = .data[[cluster_col_use]], y = score)) +
  geom_violin(scale = "width", trim = TRUE, fill = "grey85", color = "grey35") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  facet_wrap(~ signature, scales = "free_y", ncol = 3) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = "EC cluster", y = "Module score",
       title = "C3/C4 血管身份与污染程序评分")
safe_ggsave(file.path(OUTDIR, "01_identity_scores", "02_identity_module_score_violin.png"),
            p_violin, 15, 13)

# 5.2 cluster 平均分热图
cluster_score_mean <- target@meta.data %>%
  mutate(cluster = as.character(.data[[cluster_col_use]])) %>%
  group_by(cluster) %>%
  summarise(across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(cluster %in% DISPLAY_ORDER)
cluster_score_mean$cluster <- factor(cluster_score_mean$cluster, levels = DISPLAY_ORDER)
cluster_score_mean <- arrange(cluster_score_mean, cluster)
save_csv(cluster_score_mean,
         file.path(OUTDIR, "01_identity_scores", "03_cluster_mean_identity_scores.csv"))

hm_mat <- as.matrix(cluster_score_mean[, score_cols, drop = FALSE])
rownames(hm_mat) <- cluster_score_mean$cluster
hm_mat <- t(hm_mat)
png(file.path(OUTDIR, "01_identity_scores", "04_identity_module_score_heatmap.png"),
    width = 2200, height = 1800, res = 250)
pheatmap::pheatmap(
  hm_mat, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
  border_color = NA, fontsize = 10,
  main = "Identity module scores (row-scaled)"
)
dev.off()

# 5.3 marker DotPlot
marker_panels <- list(
  PanEC = signatures$Pan_endothelial,
  Arterial = signatures$Arterial_core,
  Capillary = signatures$Capillary_core,
  Angiogenic = signatures$Angiogenic_activation,
  Venous = signatures$Venous_core,
  VenousRepair = signatures$Venous_repair_like,
  ECM_EndMT = unique(c(signatures$ECM_remodeling, signatures$EndMT_regulatory)),
  Fibro_SMC = unique(c(signatures$Fibroblast, signatures$SMC_pericyte))
)
marker_features <- unique(unlist(marker_panels))
marker_features <- marker_features[marker_features %in% rownames(target)]
marker_annot <- bind_rows(lapply(names(marker_panels), function(panel) {
  data.frame(panel = panel, gene = marker_panels[[panel]], stringsAsFactors = FALSE)
})) %>% filter(gene %in% marker_features)
save_csv(marker_annot,
         file.path(OUTDIR, "01_identity_scores", "05_marker_panel_used.csv"))

p_dot <- DotPlot(target, features = marker_features, group.by = cluster_col_use,
                 assay = "RNA", dot.scale = 5) +
  RotatedAxis() +
  theme_classic(base_size = 10) +
  labs(title = "C3/C4 血管身份与污染 marker", x = NULL, y = "EC cluster")
safe_ggsave(file.path(OUTDIR, "01_identity_scores", "06_C3_C4_identity_marker_dotplot.png"),
            p_dot, 24, 6.5)

# 5.4 C3/C4 连续状态散点图
meta_sc <- target@meta.data %>% tibble::rownames_to_column("cell")
needed_sc <- c("Arterial_core_ModuleScore", "Capillary_core_ModuleScore", "Angiogenic_activation_ModuleScore",
               "Pan_endothelial_ModuleScore", "Venous_core_ModuleScore", "Venous_repair_like_ModuleScore",
               "ECM_remodeling_ModuleScore", "EndMT_regulatory_ModuleScore", "Fibroblast_ModuleScore",
               "SMC_pericyte_ModuleScore")
needed_sc <- intersect(needed_sc, colnames(meta_sc))

c3df <- meta_sc %>% filter(.data[[cluster_col_use]] == "3")
c4df <- meta_sc %>% filter(.data[[cluster_col_use]] == "4")

p_c3a <- ggplot(c3df, aes(x = Arterial_core_ModuleScore, y = Angiogenic_activation_ModuleScore)) +
  geom_point(alpha = 0.45, size = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic(base_size = 13) + labs(title = "C3：动脉程序 vs 血管生成激活")
p_c3b <- ggplot(c3df, aes(x = Arterial_core_ModuleScore, y = Capillary_core_ModuleScore)) +
  geom_point(alpha = 0.45, size = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic(base_size = 13) + labs(title = "C3：动脉程序 vs 毛细血管核心身份")
p_c4a <- ggplot(c4df, aes(x = Pan_endothelial_ModuleScore, y = ECM_remodeling_ModuleScore)) +
  geom_point(alpha = 0.45, size = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic(base_size = 13) + labs(title = "C4：内皮保留 vs ECM 重塑")
p_c4b <- ggplot(c4df, aes(x = Venous_core_ModuleScore, y = ECM_remodeling_ModuleScore)) +
  geom_point(alpha = 0.45, size = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic(base_size = 13) + labs(title = "C4：静脉核心身份 vs ECM 重塑")
p_c4c <- ggplot(c4df, aes(x = Pan_endothelial_ModuleScore, y = Fibroblast_ModuleScore)) +
  geom_point(alpha = 0.45, size = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic(base_size = 13) + labs(title = "C4：内皮保留 vs 成纤维程序")
p_c4d <- ggplot(c4df, aes(x = Arterial_core_ModuleScore, y = Venous_repair_like_ModuleScore)) +
  geom_point(alpha = 0.45, size = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_classic(base_size = 13) + labs(title = "C4：动脉程序 vs 静脉修复样程序")

safe_ggsave(file.path(OUTDIR, "01_identity_scores", "07_C3_identity_scatter.png"),
            p_c3a + p_c3b, 13, 5.8)
safe_ggsave(file.path(OUTDIR, "01_identity_scores", "08_C4_identity_scatter.png"),
            (p_c4a + p_c4b) / (p_c4c + p_c4d), 13, 11)

# 5.5 基于 C7 参考分布的启发式筛查（只用于筛查，不直接作为最终生物学结论）
ref_c7 <- meta_sc %>% filter(.data[[cluster_col_use]] == "7")
qfun <- function(x, p) as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
thresholds <- data.frame(
  metric = c("Arterial_core_ModuleScore", "Pan_endothelial_ModuleScore", "Capillary_core_ModuleScore",
             "Angiogenic_activation_ModuleScore", "Venous_core_ModuleScore", "Venous_repair_like_ModuleScore",
             "ECM_remodeling_ModuleScore", "Fibroblast_ModuleScore", "SMC_pericyte_ModuleScore"),
  q25_C7 = NA_real_, q75_C7 = NA_real_
)
for (i in seq_len(nrow(thresholds))) {
  nm <- thresholds$metric[i]
  thresholds$q25_C7[i] <- qfun(ref_c7[[nm]], 0.25)
  thresholds$q75_C7[i] <- qfun(ref_c7[[nm]], 0.75)
}
save_csv(thresholds,
         file.path(OUTDIR, "01_identity_scores", "09_C7_reference_score_thresholds.csv"))
get_thr <- function(metric, col) thresholds[thresholds$metric == metric, col][[1]]

meta_sc <- meta_sc %>%
  mutate(
    C3_transition_like = .data[[cluster_col_use]] == "3" &
      Arterial_core_ModuleScore >= get_thr("Arterial_core_ModuleScore", "q25_C7") &
      Angiogenic_activation_ModuleScore >= get_thr("Angiogenic_activation_ModuleScore", "q75_C7"),
    C3_capillary_like = .data[[cluster_col_use]] == "3" &
      Arterial_core_ModuleScore < get_thr("Arterial_core_ModuleScore", "q25_C7") &
      Capillary_core_ModuleScore >= get_thr("Capillary_core_ModuleScore", "q75_C7"),
    C4_partial_EndMT_like = .data[[cluster_col_use]] == "4" &
      Pan_endothelial_ModuleScore >= get_thr("Pan_endothelial_ModuleScore", "q25_C7") &
      ECM_remodeling_ModuleScore >= get_thr("ECM_remodeling_ModuleScore", "q75_C7"),
    C4_venous_like = .data[[cluster_col_use]] == "4" &
      Venous_core_ModuleScore >= get_thr("Venous_core_ModuleScore", "q75_C7") &
      Arterial_core_ModuleScore < get_thr("Arterial_core_ModuleScore", "q25_C7"),
    C4_contaminant_like = .data[[cluster_col_use]] == "4" &
      Pan_endothelial_ModuleScore < get_thr("Pan_endothelial_ModuleScore", "q25_C7") &
      (Fibroblast_ModuleScore >= get_thr("Fibroblast_ModuleScore", "q75_C7") |
         SMC_pericyte_ModuleScore >= get_thr("SMC_pericyte_ModuleScore", "q75_C7"))
  )

heuristic_summary <- bind_rows(
  data.frame(group = "C3", category = "transition_like",
             n = sum(meta_sc$C3_transition_like), denominator = sum(meta_sc[[cluster_col_use]] == "3")),
  data.frame(group = "C3", category = "capillary_like",
             n = sum(meta_sc$C3_capillary_like), denominator = sum(meta_sc[[cluster_col_use]] == "3")),
  data.frame(group = "C4", category = "partial_EndMT_like",
             n = sum(meta_sc$C4_partial_EndMT_like), denominator = sum(meta_sc[[cluster_col_use]] == "4")),
  data.frame(group = "C4", category = "venous_like",
             n = sum(meta_sc$C4_venous_like), denominator = sum(meta_sc[[cluster_col_use]] == "4")),
  data.frame(group = "C4", category = "contaminant_like",
             n = sum(meta_sc$C4_contaminant_like), denominator = sum(meta_sc[[cluster_col_use]] == "4"))
) %>% mutate(proportion = n / pmax(denominator, 1))
save_csv(heuristic_summary,
         file.path(OUTDIR, "01_identity_scores", "10_identity_heuristic_screening_summary.csv"))

# 把新评分写回 target 对象
saveRDS(target, file.path(OUTDIR, "01_identity_scores", "AS_EC_C1C2C3C4C5C7_identity_scored.rds"))

# --------------------------- 6. 单细胞共表达 --------------------------------
coexpression_pairs <- list(
  C3_SO17_PGF = c("SOX17", "PGF"),
  C3_EFNB2_ANGPT2 = c("EFNB2", "ANGPT2"),
  C3_GJA4_APLN = c("GJA4", "APLN"),
  C4_PECAM1_POSTN = c("PECAM1", "POSTN"),
  C4_CLDN5_COL3A1 = c("CLDN5", "COL3A1"),
  C4_VWF_SPARC = c("VWF", "SPARC"),
  C4_CDH5_COL6A2 = c("CDH5", "COL6A2"),
  C4_SOX17_POSTN = c("SOX17", "POSTN")
)
coexpression_pairs <- coexpression_pairs[vapply(coexpression_pairs, function(x) all(x %in% rownames(target)), logical(1))]
expr_detect <- get_assay_matrix(target, assay = "RNA", layer = "data")
coexpr_res <- bind_rows(lapply(names(coexpression_pairs), function(nm) {
  genes <- coexpression_pairs[[nm]]
  cl_to_use <- if (grepl("^C3_", nm)) "3" else "4"
  cc <- Cells(target)[as.character(target[[cluster_col_use]][, 1]) == cl_to_use]
  g1 <- as.numeric(expr_detect[genes[1], cc] > 0)
  g2 <- as.numeric(expr_detect[genes[2], cc] > 0)
  data.frame(
    pair = nm, cluster = cl_to_use, gene1 = genes[1], gene2 = genes[2],
    n_cells = length(cc), pct_gene1 = mean(g1) * 100, pct_gene2 = mean(g2) * 100,
    pct_double_positive = mean(g1 == 1 & g2 == 1) * 100,
    pct_double_negative = mean(g1 == 0 & g2 == 0) * 100,
    stringsAsFactors = FALSE
  )
}))
save_csv(coexpr_res,
         file.path(OUTDIR, "02_marker_coexpression", "01_C3_C4_marker_coexpression_percent.csv"))

if (nrow(coexpr_res) > 0) {
  p_co <- ggplot(coexpr_res, aes(x = reorder(pair, pct_double_positive),
                                y = pct_double_positive, fill = cluster)) +
    geom_col(width = 0.75) + coord_flip() + theme_classic(base_size = 12) +
    labs(x = NULL, y = "双阳性细胞比例（%）",
         title = "C3/C4 关键 marker 单细胞共表达")
  safe_ggsave(file.path(OUTDIR, "02_marker_coexpression", "02_marker_double_positive_barplot.png"),
              p_co, 9, 6)
}

# ---------------------- 7. 数据集/样本/组织组成检查 ---------------------------
comp_meta <- target@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(cluster = as.character(.data[[cluster_col_use]]))

# cluster 总表
cluster_n <- comp_meta %>% count(cluster, name = "n_cells")
if (!is.na(dataset_col)) {
  by_dataset <- comp_meta %>%
    count(cluster, dataset = .data[[dataset_col]], name = "n_cells") %>%
    group_by(cluster) %>% mutate(cluster_fraction = n_cells / sum(n_cells)) %>% ungroup()
  save_csv(by_dataset,
           file.path(OUTDIR, "03_composition", "01_cluster_by_dataset_counts_fractions.csv"))

  dataset_total <- comp_meta %>% count(dataset = .data[[dataset_col]], name = "dataset_total")
  cl_total <- comp_meta %>% count(cluster, name = "cluster_total")
  oe <- by_dataset %>%
    left_join(dataset_total, by = "dataset") %>%
    left_join(cl_total, by = "cluster") %>%
    mutate(
      expected = cluster_total * dataset_total / nrow(comp_meta),
      observed_expected = n_cells / pmax(expected, 1e-12),
      pearson_residual = (n_cells - expected) / sqrt(pmax(expected, 1e-12))
    )
  save_csv(oe,
           file.path(OUTDIR, "03_composition", "02_cluster_dataset_observed_expected.csv"))

  p_ds <- ggplot(by_dataset, aes(x = cluster, y = cluster_fraction, fill = as.factor(dataset))) +
    geom_col(position = "fill") + theme_classic(base_size = 12) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "EC cluster", y = "细胞组成", fill = "dataset",
         title = "各 EC cluster 的数据集组成")
  safe_ggsave(file.path(OUTDIR, "03_composition", "03_cluster_by_dataset_fraction_barplot.png"),
              p_ds, 9, 6)
}

if (!is.na(sample_col)) {
  by_sample <- comp_meta %>%
    count(cluster, sample = .data[[sample_col]], name = "n_cells") %>%
    group_by(cluster) %>% mutate(cluster_fraction = n_cells / sum(n_cells)) %>% ungroup()
  save_csv(by_sample,
           file.path(OUTDIR, "03_composition", "04_cluster_by_sample_counts_fractions.csv"))
}

if (!is.na(tissue_col)) {
  by_tissue <- comp_meta %>%
    count(cluster, tissue = .data[[tissue_col]], name = "n_cells") %>%
    group_by(cluster) %>% mutate(cluster_fraction = n_cells / sum(n_cells)) %>% ungroup()
  save_csv(by_tissue,
           file.path(OUTDIR, "03_composition", "05_cluster_by_tissue_counts_fractions.csv"))
}

if (!is.na(disease_col)) {
  by_disease <- comp_meta %>%
    count(cluster, disease = .data[[disease_col]], name = "n_cells") %>%
    group_by(cluster) %>% mutate(cluster_fraction = n_cells / sum(n_cells)) %>% ungroup()
  save_csv(by_disease,
           file.path(OUTDIR, "03_composition", "06_cluster_by_disease_counts_fractions.csv"))
}

cluster_support <- cluster_n
if (!is.na(dataset_col)) {
  tmp <- comp_meta %>% group_by(cluster) %>%
    summarise(n_datasets = n_distinct(.data[[dataset_col]]),
              largest_dataset_fraction = max(table(.data[[dataset_col]])) / n(), .groups = "drop")
  cluster_support <- left_join(cluster_support, tmp, by = "cluster")
}
if (!is.na(sample_col)) {
  tmp <- comp_meta %>% group_by(cluster) %>%
    summarise(n_samples = n_distinct(.data[[sample_col]]),
              largest_sample_fraction = max(table(.data[[sample_col]])) / n(), .groups = "drop")
  cluster_support <- left_join(cluster_support, tmp, by = "cluster")
}
save_csv(cluster_support,
         file.path(OUTDIR, "03_composition", "07_cluster_cross_dataset_sample_support.csv"))

# 样本级平均分，避免仅由大量细胞样本主导
if (!is.na(sample_col)) {
  sample_score_mean <- target@meta.data %>%
    mutate(cluster = as.character(.data[[cluster_col_use]]),
           sample = as.character(.data[[sample_col]])) %>%
    group_by(sample, cluster) %>%
    summarise(across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
              n_cells = n(), .groups = "drop")
  if (!is.na(dataset_col)) {
    sample_to_dataset <- target@meta.data %>%
      transmute(sample = as.character(.data[[sample_col]]),
                dataset = as.character(.data[[dataset_col]])) %>% distinct()
    sample_score_mean <- left_join(sample_score_mean, sample_to_dataset, by = "sample")
  }
  save_csv(sample_score_mean,
           file.path(OUTDIR, "03_composition", "08_sample_level_identity_scores.csv"))
}

# ---------------------------- 8. 局部重分析 ---------------------------------
cat("\nStarting local PCA/Harmony reanalysis...\n")
local <- prepare_local_object(target, dataset_col = dataset_col, run_harmony = TRUE)

# UMAP 图：原 EC cluster / 新局部 cluster / dataset
p_upca_orig <- DimPlot(local, reduction = "umap_pca_local", group.by = cluster_col_use,
                       label = TRUE, repel = TRUE, raster = TRUE) +
  ggtitle("RNA-PCA local UMAP：原 EC cluster")
p_upca_new <- DimPlot(local, reduction = "umap_pca_local", group.by = "local_cluster_pca",
                      label = TRUE, repel = TRUE, raster = TRUE) +
  ggtitle("RNA-PCA local UMAP：重新聚类")
if (!is.na(dataset_col)) {
  p_upca_ds <- DimPlot(local, reduction = "umap_pca_local", group.by = dataset_col,
                       raster = TRUE) + ggtitle("RNA-PCA local UMAP：dataset")
  p_upca_all <- p_upca_orig + p_upca_new + p_upca_ds + plot_layout(ncol = 3)
} else {
  p_upca_all <- p_upca_orig + p_upca_new
}
safe_ggsave(file.path(OUTDIR, "04_local_reanalysis", "01_local_RNA_PCA_UMAP.png"),
            p_upca_all, ifelse(is.na(dataset_col), 13, 19), 6.5)

has_harmony <- "harmony_local" %in% names(local@reductions)
if (has_harmony) {
  p_uh_orig <- DimPlot(local, reduction = "umap_harmony_local", group.by = cluster_col_use,
                       label = TRUE, repel = TRUE, raster = TRUE) +
    ggtitle("Harmony local UMAP：原 EC cluster")
  p_uh_new <- DimPlot(local, reduction = "umap_harmony_local", group.by = "local_cluster_harmony",
                      label = TRUE, repel = TRUE, raster = TRUE) +
    ggtitle("Harmony local UMAP：重新聚类")
  if (!is.na(dataset_col)) {
    p_uh_ds <- DimPlot(local, reduction = "umap_harmony_local", group.by = dataset_col,
                       raster = TRUE) + ggtitle("Harmony local UMAP：dataset")
    p_uh_all <- p_uh_orig + p_uh_new + p_uh_ds + plot_layout(ncol = 3)
  } else {
    p_uh_all <- p_uh_orig + p_uh_new
  }
  safe_ggsave(file.path(OUTDIR, "04_local_reanalysis", "02_local_Harmony_UMAP.png"),
              p_uh_all, ifelse(is.na(dataset_col), 13, 19), 6.5)
}

# 新局部 cluster 与原 cluster 的映射
map_pca <- local@meta.data %>%
  transmute(local_cluster = as.character(local_cluster_pca),
            original_cluster = as.character(.data[[cluster_col_use]])) %>%
  count(local_cluster, original_cluster, name = "n") %>%
  group_by(local_cluster) %>% mutate(proportion = n / sum(n)) %>% ungroup()
save_csv(map_pca,
         file.path(OUTDIR, "04_local_reanalysis", "03_local_PCA_cluster_mapping.csv"))

if (has_harmony) {
  map_h <- local@meta.data %>%
    transmute(local_cluster = as.character(local_cluster_harmony),
              original_cluster = as.character(.data[[cluster_col_use]])) %>%
    count(local_cluster, original_cluster, name = "n") %>%
    group_by(local_cluster) %>% mutate(proportion = n / sum(n)) %>% ungroup()
  save_csv(map_h,
           file.path(OUTDIR, "04_local_reanalysis", "04_local_Harmony_cluster_mapping.csv"))
}

# kNN/SNN cluster 邻接关系
adj_pca <- aggregate_graph_by_cluster(local, "pca_snn", cluster_col_use)
if (!is.null(adj_pca)) {
  save_csv(cbind(cluster = rownames(adj_pca$raw), as.data.frame(adj_pca$raw)),
           file.path(OUTDIR, "04_local_reanalysis", "05_PCA_cluster_graph_adjacency_raw.csv"))
  save_csv(cbind(cluster = rownames(adj_pca$row_prop), as.data.frame(adj_pca$row_prop)),
           file.path(OUTDIR, "04_local_reanalysis", "06_PCA_cluster_graph_adjacency_rowprop.csv"))
  png(file.path(OUTDIR, "04_local_reanalysis", "07_PCA_cluster_graph_adjacency_heatmap.png"),
      width = 1800, height = 1600, res = 250)
  pheatmap(adj_pca$row_prop, cluster_rows = FALSE, cluster_cols = FALSE,
           main = "PCA SNN adjacency (row proportion)", border_color = NA)
  dev.off()
}
if (has_harmony) {
  adj_h <- aggregate_graph_by_cluster(local, "harmony_snn", cluster_col_use)
  if (!is.null(adj_h)) {
    save_csv(cbind(cluster = rownames(adj_h$raw), as.data.frame(adj_h$raw)),
             file.path(OUTDIR, "04_local_reanalysis", "08_Harmony_cluster_graph_adjacency_raw.csv"))
    save_csv(cbind(cluster = rownames(adj_h$row_prop), as.data.frame(adj_h$row_prop)),
             file.path(OUTDIR, "04_local_reanalysis", "09_Harmony_cluster_graph_adjacency_rowprop.csv"))
    png(file.path(OUTDIR, "04_local_reanalysis", "10_Harmony_cluster_graph_adjacency_heatmap.png"),
        width = 1800, height = 1600, res = 250)
    pheatmap(adj_h$row_prop, cluster_rows = FALSE, cluster_cols = FALSE,
             main = "Harmony SNN adjacency (row proportion)", border_color = NA)
    dev.off()
  }
}

saveRDS(local,
        file.path(OUTDIR, "04_local_reanalysis", "AS_EC_C1C2C3C4C5C7_local_reanalysis.rds"))

# --------------------------- 9. Slingshot 验证 -------------------------------
cat("\nRunning Slingshot without forced terminal clusters...\n")
trajectory_all <- list()

# 9.1 原 cluster + PCA
emb_pca <- Embeddings(local, "pca")[, DIMS_USE, drop = FALSE]
labs_orig <- as.character(local[[cluster_col_use]][, 1])
sds_pca_orig <- run_slingshot_matrix(emb_pca, labs_orig, ROOT_CLUSTER)
sum_pca_orig <- lineage_summary(sds_pca_orig)
sum_pca_orig$space <- "PCA"
sum_pca_orig$cluster_definition <- "original_EC_cluster"
save_csv(sum_pca_orig,
         file.path(OUTDIR, "05_slingshot", "01_PCA_original_cluster_lineages.csv"))
plot_slingshot_2d(emb_pca, labs_orig, sds_pca_orig,
                  file.path(OUTDIR, "05_slingshot", "02_PCA_original_cluster_slingshot.png"),
                  "Slingshot on RNA-PCA, original EC clusters, root=C7")
local <- add_pseudotime_and_plot(local, sds_pca_orig, "umap_pca_local",
                                 "PCA_original", file.path(OUTDIR, "05_slingshot"))
trajectory_all[["PCA_original"]] <- sum_pca_orig
saveRDS(sds_pca_orig,
        file.path(OUTDIR, "05_slingshot", "PCA_original_cluster_slingshot.rds"))

# 9.2 原 cluster + Harmony
if (has_harmony) {
  emb_h <- Embeddings(local, "harmony_local")[, DIMS_USE, drop = FALSE]
  sds_h_orig <- run_slingshot_matrix(emb_h, labs_orig, ROOT_CLUSTER)
  sum_h_orig <- lineage_summary(sds_h_orig)
  sum_h_orig$space <- "Harmony"
  sum_h_orig$cluster_definition <- "original_EC_cluster"
  save_csv(sum_h_orig,
           file.path(OUTDIR, "05_slingshot", "03_Harmony_original_cluster_lineages.csv"))
  plot_slingshot_2d(emb_h, labs_orig, sds_h_orig,
                    file.path(OUTDIR, "05_slingshot", "04_Harmony_original_cluster_slingshot.png"),
                    "Slingshot on Harmony, original EC clusters, root=C7")
  local <- add_pseudotime_and_plot(local, sds_h_orig, "umap_harmony_local",
                                   "Harmony_original", file.path(OUTDIR, "05_slingshot"))
  trajectory_all[["Harmony_original"]] <- sum_h_orig
  saveRDS(sds_h_orig,
          file.path(OUTDIR, "05_slingshot", "Harmony_original_cluster_slingshot.rds"))
}

# 9.3 局部重新聚类 + PCA：根簇由 C7 占比最高的局部 cluster 自动确定
root_local_pca <- map_pca %>%
  filter(original_cluster == "7") %>% arrange(desc(proportion), desc(n)) %>% slice(1) %>% pull(local_cluster)
lab_local_pca <- as.character(local$local_cluster_pca)
sds_pca_local <- run_slingshot_matrix(emb_pca, lab_local_pca, root_local_pca)
sum_pca_local <- lineage_summary(sds_pca_local)
sum_pca_local$space <- "PCA"
sum_pca_local$cluster_definition <- "local_recluster"
sum_pca_local$root_local_cluster <- root_local_pca
save_csv(sum_pca_local,
         file.path(OUTDIR, "05_slingshot", "05_PCA_local_recluster_lineages.csv"))
plot_slingshot_2d(emb_pca, lab_local_pca, sds_pca_local,
                  file.path(OUTDIR, "05_slingshot", "06_PCA_local_recluster_slingshot.png"),
                  paste0("Slingshot on RNA-PCA local clusters, root local cluster=", root_local_pca))
saveRDS(sds_pca_local,
        file.path(OUTDIR, "05_slingshot", "PCA_local_recluster_slingshot.rds"))

# 将局部谱系按局部 cluster 的多数原标签进行映射
major_map_pca <- map_pca %>% group_by(local_cluster) %>% slice_max(proportion, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(local_cluster, majority_original = original_cluster)
map_local_lineages <- function(summary_df, major_map) {
  if (nrow(summary_df) == 0) return(summary_df)
  summary_df %>% rowwise() %>% mutate(
    mapped_original_sequence = {
      seq_local <- strsplit(sequence, "->", fixed = TRUE)[[1]]
      mapped <- major_map$majority_original[match(seq_local, major_map$local_cluster)]
      paste(mapped, collapse = "->")
    }
  ) %>% ungroup()
}
mapped_pca_local <- map_local_lineages(sum_pca_local, major_map_pca)
save_csv(mapped_pca_local,
         file.path(OUTDIR, "05_slingshot", "07_PCA_local_lineages_mapped_to_original.csv"))

if (has_harmony) {
  root_local_h <- map_h %>%
    filter(original_cluster == "7") %>% arrange(desc(proportion), desc(n)) %>% slice(1) %>% pull(local_cluster)
  lab_local_h <- as.character(local$local_cluster_harmony)
  sds_h_local <- run_slingshot_matrix(emb_h, lab_local_h, root_local_h)
  sum_h_local <- lineage_summary(sds_h_local)
  sum_h_local$space <- "Harmony"
  sum_h_local$cluster_definition <- "local_recluster"
  sum_h_local$root_local_cluster <- root_local_h
  save_csv(sum_h_local,
           file.path(OUTDIR, "05_slingshot", "08_Harmony_local_recluster_lineages.csv"))
  plot_slingshot_2d(emb_h, lab_local_h, sds_h_local,
                    file.path(OUTDIR, "05_slingshot", "09_Harmony_local_recluster_slingshot.png"),
                    paste0("Slingshot on Harmony local clusters, root local cluster=", root_local_h))
  saveRDS(sds_h_local,
          file.path(OUTDIR, "05_slingshot", "Harmony_local_recluster_slingshot.rds"))
  major_map_h <- map_h %>% group_by(local_cluster) %>% slice_max(proportion, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(local_cluster, majority_original = original_cluster)
  mapped_h_local <- map_local_lineages(sum_h_local, major_map_h)
  save_csv(mapped_h_local,
           file.path(OUTDIR, "05_slingshot", "10_Harmony_local_lineages_mapped_to_original.csv"))
}

full_traj_summary <- bind_rows(trajectory_all)
save_csv(full_traj_summary,
         file.path(OUTDIR, "05_slingshot", "11_full_trajectory_summary.csv"))

full_flags <- bind_rows(lapply(names(trajectory_all), function(nm) {
  cbind(model = nm, trajectory_flags(trajectory_all[[nm]]))
}))
save_csv(full_flags,
         file.path(OUTDIR, "05_slingshot", "12_full_trajectory_recovery_flags.csv"))

saveRDS(local,
        file.path(OUTDIR, "05_slingshot", "AS_EC_local_reanalysis_with_new_pseudotime.rds"))

# ---------------------- 10. leave-one-dataset-out ----------------------------
lodo_results <- data.frame()
if (RUN_LODO && !is.na(dataset_col) && length(unique(local[[dataset_col]][, 1])) >= 3) {
  datasets <- sort(unique(as.character(local[[dataset_col]][, 1])))
  for (drop_ds in datasets) {
    cat("LODO: dropping ", drop_ds, "\n", sep = "")
    keep_cells <- Cells(target)[as.character(target[[dataset_col]][, 1]) != drop_ds]
    sub <- subset(target, cells = keep_cells)
    row <- tryCatch({
      sub <- prepare_local_object(sub, dataset_col = dataset_col, run_harmony = TRUE)
      labs <- as.character(sub[[cluster_col_use]][, 1])
      pca_emb <- Embeddings(sub, "pca")[, DIMS_USE, drop = FALSE]
      sds_p <- run_slingshot_matrix(pca_emb, labs, ROOT_CLUSTER)
      sum_p <- lineage_summary(sds_p)
      flag_p <- trajectory_flags(sum_p)
      line_p <- paste(sum_p$sequence, collapse = " | ")

      if ("harmony_local" %in% names(sub@reductions)) {
        h_emb <- Embeddings(sub, "harmony_local")[, DIMS_USE, drop = FALSE]
        sds_h <- run_slingshot_matrix(h_emb, labs, ROOT_CLUSTER)
        sum_h <- lineage_summary(sds_h)
        flag_h <- trajectory_flags(sum_h)
        line_h <- paste(sum_h$sequence, collapse = " | ")
      } else {
        flag_h <- data.frame(n_lineages = NA, L1_exact = NA, L1_core = NA,
                             L2 = NA, both_L1core_L2 = NA)
        line_h <- NA_character_
      }

      data.frame(
        omitted_dataset = drop_ds,
        n_cells = ncol(sub),
        n_samples = if (!is.na(sample_col)) length(unique(sub[[sample_col]][, 1])) else NA_integer_,
        PCA_n_lineages = flag_p$n_lineages,
        PCA_L1_exact = flag_p$L1_exact,
        PCA_L1_core = flag_p$L1_core,
        PCA_L2 = flag_p$L2,
        PCA_both = flag_p$both_L1core_L2,
        PCA_sequences = line_p,
        Harmony_n_lineages = flag_h$n_lineages,
        Harmony_L1_exact = flag_h$L1_exact,
        Harmony_L1_core = flag_h$L1_core,
        Harmony_L2 = flag_h$L2,
        Harmony_both = flag_h$both_L1core_L2,
        Harmony_sequences = line_h,
        error = NA_character_,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        omitted_dataset = drop_ds, n_cells = length(keep_cells), n_samples = NA_integer_,
        PCA_n_lineages = NA, PCA_L1_exact = NA, PCA_L1_core = NA, PCA_L2 = NA,
        PCA_both = NA, PCA_sequences = NA_character_, Harmony_n_lineages = NA,
        Harmony_L1_exact = NA, Harmony_L1_core = NA, Harmony_L2 = NA,
        Harmony_both = NA, Harmony_sequences = NA_character_,
        error = conditionMessage(e), stringsAsFactors = FALSE
      )
    })
    lodo_results <- bind_rows(lodo_results, row)
    save_csv(lodo_results,
             file.path(OUTDIR, "06_LODO", "01_leave_one_dataset_out_branch_recovery.csv"))
    gc()
  }

  lodo_summary <- data.frame(
    metric = c("PCA_L1_core", "PCA_L2", "PCA_both",
               "Harmony_L1_core", "Harmony_L2", "Harmony_both"),
    recovery_rate = c(
      mean(lodo_results$PCA_L1_core, na.rm = TRUE),
      mean(lodo_results$PCA_L2, na.rm = TRUE),
      mean(lodo_results$PCA_both, na.rm = TRUE),
      mean(lodo_results$Harmony_L1_core, na.rm = TRUE),
      mean(lodo_results$Harmony_L2, na.rm = TRUE),
      mean(lodo_results$Harmony_both, na.rm = TRUE)
    )
  )
  save_csv(lodo_summary,
           file.path(OUTDIR, "06_LODO", "02_LODO_recovery_rate_summary.csv"))
}

# ---------------------- 11. 样本级 bootstrap ---------------------------------
bootstrap_results <- data.frame()
if (RUN_BOOTSTRAP && !is.na(sample_col) && length(unique(local[[sample_col]][, 1])) >= 10) {
  samples_all <- sort(unique(as.character(local[[sample_col]][, 1])))
  n_take <- max(8, floor(length(samples_all) * BOOT_FRAC))
  emb_pca_full <- Embeddings(local, "pca")[, DIMS_USE, drop = FALSE]
  emb_h_full <- if (has_harmony) Embeddings(local, "harmony_local")[, DIMS_USE, drop = FALSE] else NULL

  for (b in seq_len(N_BOOT)) {
    set.seed(20260728 + b)
    selected_samples <- sample(samples_all, size = n_take, replace = FALSE)
    keep <- as.character(local[[sample_col]][, 1]) %in% selected_samples
    labs_b <- as.character(local[[cluster_col_use]][keep, 1])

    row <- tryCatch({
      if (!all(TARGET_CLUSTERS %in% unique(labs_b))) {
        stop("Not all target clusters represented in this bootstrap.")
      }
      sds_p <- run_slingshot_matrix(emb_pca_full[keep, , drop = FALSE], labs_b, ROOT_CLUSTER)
      sum_p <- lineage_summary(sds_p)
      flag_p <- trajectory_flags(sum_p)
      line_p <- paste(sum_p$sequence, collapse = " | ")

      if (!is.null(emb_h_full)) {
        sds_h <- run_slingshot_matrix(emb_h_full[keep, , drop = FALSE], labs_b, ROOT_CLUSTER)
        sum_h <- lineage_summary(sds_h)
        flag_h <- trajectory_flags(sum_h)
        line_h <- paste(sum_h$sequence, collapse = " | ")
      } else {
        flag_h <- data.frame(n_lineages = NA, L1_exact = NA, L1_core = NA,
                             L2 = NA, both_L1core_L2 = NA)
        line_h <- NA_character_
      }

      data.frame(
        bootstrap = b, n_samples = length(selected_samples), n_cells = sum(keep),
        PCA_L1_exact = flag_p$L1_exact, PCA_L1_core = flag_p$L1_core,
        PCA_L2 = flag_p$L2, PCA_both = flag_p$both_L1core_L2,
        PCA_sequences = line_p,
        Harmony_L1_exact = flag_h$L1_exact, Harmony_L1_core = flag_h$L1_core,
        Harmony_L2 = flag_h$L2, Harmony_both = flag_h$both_L1core_L2,
        Harmony_sequences = line_h, error = NA_character_, stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        bootstrap = b, n_samples = length(selected_samples), n_cells = sum(keep),
        PCA_L1_exact = NA, PCA_L1_core = NA, PCA_L2 = NA, PCA_both = NA,
        PCA_sequences = NA_character_, Harmony_L1_exact = NA,
        Harmony_L1_core = NA, Harmony_L2 = NA, Harmony_both = NA,
        Harmony_sequences = NA_character_, error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
    bootstrap_results <- bind_rows(bootstrap_results, row)
    if (b %% 5 == 0 || b == N_BOOT) {
      cat("Bootstrap ", b, "/", N_BOOT, " complete\n", sep = "")
      save_csv(bootstrap_results,
               file.path(OUTDIR, "07_bootstrap", "01_sample_bootstrap_branch_recovery.csv"))
    }
  }

  bootstrap_summary <- data.frame(
    metric = c("PCA_L1_core", "PCA_L2", "PCA_both",
               "Harmony_L1_core", "Harmony_L2", "Harmony_both"),
    recovery_rate = c(
      mean(bootstrap_results$PCA_L1_core, na.rm = TRUE),
      mean(bootstrap_results$PCA_L2, na.rm = TRUE),
      mean(bootstrap_results$PCA_both, na.rm = TRUE),
      mean(bootstrap_results$Harmony_L1_core, na.rm = TRUE),
      mean(bootstrap_results$Harmony_L2, na.rm = TRUE),
      mean(bootstrap_results$Harmony_both, na.rm = TRUE)
    ),
    n_success = c(
      sum(!is.na(bootstrap_results$PCA_L1_core)),
      sum(!is.na(bootstrap_results$PCA_L2)),
      sum(!is.na(bootstrap_results$PCA_both)),
      sum(!is.na(bootstrap_results$Harmony_L1_core)),
      sum(!is.na(bootstrap_results$Harmony_L2)),
      sum(!is.na(bootstrap_results$Harmony_both))
    )
  )
  save_csv(bootstrap_summary,
           file.path(OUTDIR, "07_bootstrap", "02_bootstrap_recovery_rate_summary.csv"))
}

# --------------------------- 12. 最终汇总指标 --------------------------------
# 这里输出“证据指标”，不自动替你下最终生物学结论。
identity_metrics <- heuristic_summary %>%
  select(group, category, n, denominator, proportion)

support_C3 <- cluster_support %>% filter(cluster == "3")
support_C4 <- cluster_support %>% filter(cluster == "4")

final_metrics <- list(
  identity_heuristics = identity_metrics,
  cross_dataset_support_C3 = support_C3,
  cross_dataset_support_C4 = support_C4,
  full_trajectory_flags = full_flags,
  LODO = if (exists("lodo_summary")) lodo_summary else data.frame(),
  bootstrap = if (exists("bootstrap_summary")) bootstrap_summary else data.frame()
)

saveRDS(final_metrics,
        file.path(OUTDIR, "08_final_summary", "01_all_decision_metrics.rds"))

# 汇总成一个长表
final_long <- bind_rows(
  identity_metrics %>% transmute(section = "identity", metric = paste(group, category, sep = "_"),
                                 value = proportion, detail = paste0(n, "/", denominator)),
  full_flags %>% pivot_longer(cols = -model, names_to = "metric", values_to = "value") %>%
    mutate(section = "full_trajectory", detail = model) %>%
    select(section, metric, value, detail),
  if (exists("lodo_summary")) lodo_summary %>%
    transmute(section = "LODO", metric, value = recovery_rate, detail = NA_character_) else NULL,
  if (exists("bootstrap_summary")) bootstrap_summary %>%
    transmute(section = "bootstrap", metric, value = recovery_rate,
              detail = paste0("n_success=", n_success)) else NULL
)
save_csv(final_long,
         file.path(OUTDIR, "08_final_summary", "02_final_evidence_metrics.csv"))

interpretation_text <- c(
  "本脚本不会仅凭一个自动阈值直接宣布双分支成立或不成立。请重点检查以下文件：",
  "1. 01_identity_scores/04_identity_module_score_heatmap.png",
  "2. 01_identity_scores/07_C3_identity_scatter.png",
  "3. 01_identity_scores/08_C4_identity_scatter.png",
  "4. 03_composition/07_cluster_cross_dataset_sample_support.csv",
  "5. 04_local_reanalysis/01_local_RNA_PCA_UMAP.png",
  "6. 04_local_reanalysis/02_local_Harmony_UMAP.png",
  "7. 05_slingshot/12_full_trajectory_recovery_flags.csv",
  "8. 06_LODO/02_LODO_recovery_rate_summary.csv",
  "9. 07_bootstrap/02_bootstrap_recovery_rate_summary.csv",
  "",
  "主要判断原则：",
  "- C3 若保留动脉程序，同时血管生成激活上升，但毛细血管核心身份并未完全占主导，才支持共同过渡态。",
  "- C4 若仍保留 Pan-EC，同时 ECM/EndMT 上升，且静脉核心、成纤维/周细胞污染不能解释其主要信号，才支持 partial EndMT-like。",
  "- 两条分支应在 PCA 与 Harmony 中均能自然恢复，并在删除单个数据集及样本级 bootstrap 中保持较高恢复率。",
  "- 如果短分支只在 Harmony、单一数据集或少数样本中出现，应视为高风险，而不是直接作为论文核心结论。"
)
writeLines(interpretation_text,
           file.path(OUTDIR, "08_final_summary", "03_HOW_TO_INTERPRET_RESULTS.txt"))

cat("\nAll analyses completed.\n")
cat("End time:", format(Sys.time()), "\n")
cat("Main output directory:", OUTDIR, "\n")
