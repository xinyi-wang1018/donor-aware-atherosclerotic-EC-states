#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# GSE159677 独立外部验证：严格版一键重跑脚本
#
# 核心修正：
# 1. 修复 Harmony 参数：明确使用 reduction.use = "pca"；按 patient_id 校正个体差异，Harmony 失败即停止。
# 2. 全图谱先做 cluster 级纯度审计，再做 cell 级严格 EC 过滤。
# 3. 排除明显免疫、平滑肌/周细胞、成纤维细胞污染。
# 4. 将 C7 拆为“动脉身份”和“保护性稳态”；将 C3 固定为 ANGPT2/PGF 激活核心。
# 5. C4 使用不含 COL1A1/COL1A2/COL3A1 的主评分，并增加广义 ECM 敏感性分析。
# 6. 所有 AC-PA 比较均按 3 位患者配对；同时输出 pseudobulk edgeR。
# 7. 自动生成小型审阅包，不再要求上传大型 RDS。
###############################################################################

OLD_BASE="/root/GSE159677_external_validation"
BASE="/root/GSE159677_external_validation_strict"
INPUT_DIR="${BASE}/input"
AGG_DIR="${BASE}/aggregate_matrix"
OUT_DIR="${BASE}/results"
OBJ_DIR="${BASE}/objects"
LOG_DIR="${BASE}/logs"
REF_EC="/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
R_SCRIPT="${BASE}/GSE159677_external_validation_STRICT.R"

mkdir -p "$INPUT_DIR" "$LOG_DIR"
# 每次重跑清理旧的中间矩阵、结果和对象，避免旧文件混入；保留 input 和历史日志。
rm -rf "$AGG_DIR" "$OUT_DIR" "$OBJ_DIR"
mkdir -p "$AGG_DIR" "$OUT_DIR" "$OBJ_DIR"

exec > >(tee -a "${LOG_DIR}/run_$(date +%Y%m%d_%H%M%S).log") 2>&1

echo "============================================================"
echo "GSE159677 strict independent external validation"
echo "Start: $(date)"
echo "BASE: $BASE"
echo "============================================================"

# -----------------------------------------------------------------------------
# 1. 自动寻找你已经上传的 293.5 MB 聚合矩阵；不需要重新上传
# -----------------------------------------------------------------------------
mapfile -t AGG_ARCHIVES < <(
  find "$INPUT_DIR" "$OLD_BASE/input" -maxdepth 2 -type f 2>/dev/null \
    \( -iname 'GSE159677_AGGREGATEMAPPED-tisCAR6samples_featurebcmatrixfiltered.tar.gz' \
       -o -iname 'GSE159677_AGGREGATEMAPPED-tisCAR6samples_featurebcmatrixfiltered.tgz' \
       -o -iname '*AGGREGATEMAPPED*tisCAR6samples*featurebcmatrixfiltered*.tar.gz' \
       -o -iname '*AGGREGATEMAPPED*tisCAR6samples*featurebcmatrixfiltered*.tgz' \) \
    | sort -u
)

if (( ${#AGG_ARCHIVES[@]} < 1 )); then
  echo "ERROR: 没有找到 GSE159677 聚合矩阵。" >&2
  echo "应存在以下文件：" >&2
  echo "${OLD_BASE}/input/GSE159677_AGGREGATEMAPPED-tisCAR6samples_featurebcmatrixfiltered.tar.gz" >&2
  exit 1
fi

# 如果旧目录和新目录各有一份，优先使用旧目录现有文件，避免重复上传。
AGG_ARCHIVE="${AGG_ARCHIVES[0]}"
for candidate in "${AGG_ARCHIVES[@]}"; do
  if [[ "$candidate" == "$OLD_BASE/input/"* ]]; then
    AGG_ARCHIVE="$candidate"
    break
  fi
done

echo "Using aggregate archive: $AGG_ARCHIVE"

# -----------------------------------------------------------------------------
# 2. 解压数据
# -----------------------------------------------------------------------------
tar -xzf "$AGG_ARCHIVE" -C "$AGG_DIR"

n_matrix=$(find "$AGG_DIR" -type f -iname '*matrix.mtx*' | wc -l)
n_feature=$(find "$AGG_DIR" -type f \( -iname '*features.tsv*' -o -iname '*genes.tsv*' \) | wc -l)
n_barcode=$(find "$AGG_DIR" -type f -iname '*barcodes.tsv*' | wc -l)
if [[ "$n_matrix" -lt 1 || "$n_feature" -lt 1 || "$n_barcode" -lt 1 ]]; then
  echo "ERROR: 解压后没有找到 matrix/features/barcodes 三件套。" >&2
  find "$AGG_DIR" -type f -printf '%p\n' | sort
  exit 1
fi

# -----------------------------------------------------------------------------
# 3. 固定 6 个 GEM 后缀对应的 3 对患者
# -----------------------------------------------------------------------------
cat > "${BASE}/sample_metadata.csv" <<'CSV'
gem_suffix,sample_id,patient_id,condition,tissue,GSM
1,P1_PA,P1,PA,proximal_adjacent,GSM4837524
2,P1_AC,P1,AC,atherosclerotic_core,GSM4837523
3,P2_PA,P2,PA,proximal_adjacent,GSM4837526
4,P2_AC,P2,AC,atherosclerotic_core,GSM4837525
5,P3_PA,P3,PA,proximal_adjacent,GSM4837528
6,P3_AC,P3,AC,atherosclerotic_core,GSM4837527
CSV

# -----------------------------------------------------------------------------
# 4. 写入严格版 R 分析脚本
# -----------------------------------------------------------------------------
cat > "$R_SCRIPT" <<'RSCRIPT'
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 40 * 1024^3)
set.seed(20260731)

BASE <- "/root/GSE159677_external_validation_strict"
AGG_DIR <- file.path(BASE, "aggregate_matrix")
OUT_DIR <- file.path(BASE, "results")
OBJ_DIR <- file.path(BASE, "objects")
REF_EC <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr", "patchwork", "pheatmap", "harmony")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required R packages: ", paste(missing, collapse = ", "),
    "\n请进入之前可运行 Seurat/Harmony 的 R 环境后重跑。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
})

white_theme <- theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )
theme_set(white_theme)

save_png <- function(plot, filename, width = 10, height = 7, dpi = 300) {
  ggsave(
    filename = file.path(OUT_DIR, filename), plot = plot,
    width = width, height = height, dpi = dpi,
    bg = "white", limitsize = FALSE
  )
}

join_layers_safe <- function(obj, assay = "RNA") {
  obj2 <- tryCatch(
    SeuratObject::JoinLayers(obj, assay = assay),
    error = function(e1) tryCatch(
      Seurat::JoinLayers(obj, assay = assay),
      error = function(e2) obj
    )
  )
  obj2
}

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

open_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else if (grepl("\\.(bz2|bz)$", path, ignore.case = TRUE)) {
    bzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

find_one_file <- function(root, patterns) {
  ff <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  for (pat in patterns) {
    hit <- ff[grepl(pat, basename(ff), ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  stop("Could not find required file under ", root, ": ", paste(patterns, collapse = " | "))
}

read_10x_any_compression <- function(sample_root) {
  matrix_file <- find_one_file(sample_root, c("matrix\\.mtx(\\.gz|\\.bz2|\\.bz)?$"))
  feature_file <- find_one_file(sample_root, c(
    "features\\.tsv(\\.gz|\\.bz2|\\.bz)?$",
    "genes\\.tsv(\\.gz|\\.bz2|\\.bz)?$"
  ))
  barcode_file <- find_one_file(sample_root, c("barcodes\\.tsv(\\.gz|\\.bz2|\\.bz)?$"))

  message("Reading aggregate matrix:")
  message("  matrix:   ", matrix_file)
  message("  features: ", feature_file)
  message("  barcodes: ", barcode_file)

  cm <- open_text_connection(matrix_file)
  mat <- Matrix::readMM(cm)
  close(cm)

  cf <- open_text_connection(feature_file)
  features <- utils::read.delim(
    cf, header = FALSE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE
  )
  close(cf)

  cb <- open_text_connection(barcode_file)
  barcodes <- utils::read.delim(
    cb, header = FALSE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE
  )
  close(cb)

  if (nrow(features) != nrow(mat)) stop("Feature rows do not equal matrix rows.")
  if (nrow(barcodes) != ncol(mat)) stop("Barcode rows do not equal matrix columns.")

  gene_col <- if (ncol(features) >= 2) 2 else 1
  gene_names <- as.character(features[[gene_col]])
  bad_gene <- is.na(gene_names) | gene_names == ""
  gene_names[bad_gene] <- as.character(features[[1]])[bad_gene]
  gene_names <- make.unique(gene_names)

  mat <- as(mat, "CsparseMatrix")
  rownames(mat) <- gene_names
  colnames(mat) <- as.character(barcodes[[1]])
  mat
}

run_harmony_strict <- function(obj, group_var, reduction_name = "pca", save_name = "harmony") {
  message("Running Harmony by ", group_var, " ...")
  out <- tryCatch(
    harmony::RunHarmony(
      object = obj,
      group.by.vars = group_var,
      reduction.use = reduction_name,
      reduction.save = save_name,
      verbose = TRUE
    ),
    error = function(e1) {
      message("First Harmony call failed: ", conditionMessage(e1))
      tryCatch(
        harmony::RunHarmony(
          object = obj,
          group.by.vars = group_var,
          reduction.use = reduction_name,
          verbose = TRUE
        ),
        error = function(e2) {
          stop(
            "Harmony failed and strict analysis will not fall back to PCA.\n",
            "First error: ", conditionMessage(e1), "\n",
            "Second error: ", conditionMessage(e2)
          )
        }
      )
    }
  )
  if (!save_name %in% Reductions(out) && !"harmony" %in% Reductions(out)) {
    stop("Harmony returned without a harmony reduction.")
  }
  out
}

marker_count <- function(count_mat, genes) {
  genes <- intersect(genes, rownames(count_mat))
  if (length(genes) == 0) return(rep(0L, ncol(count_mat)))
  as.integer(Matrix::colSums(count_mat[genes, , drop = FALSE] > 0))
}

mean_log_score <- function(data_mat, genes) {
  genes <- intersect(genes, rownames(data_mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(data_mat)))
  as.numeric(Matrix::colMeans(data_mat[genes, , drop = FALSE]))
}

gene_z_score <- function(data_mat, genes, cap = 5) {
  genes <- intersect(genes, rownames(data_mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(data_mat)))
  x <- as.matrix(data_mat[genes, , drop = FALSE])
  mu <- rowMeans(x)
  s <- apply(x, 1, stats::sd)
  s[!is.finite(s) | s < 1e-8] <- 1
  z <- sweep(x, 1, mu, "-")
  z <- sweep(z, 1, s, "/")
  z[z > cap] <- cap
  z[z < -cap] <- -cap
  colMeans(z)
}

paired_stats <- function(df, value_col, feature_col, expected_map = NULL) {
  features <- unique(df[[feature_col]])
  out <- lapply(features, function(feat) {
    dd <- df[df[[feature_col]] == feat, c("patient_id", "condition", value_col), drop = FALSE]
    wide <- tidyr::pivot_wider(dd, names_from = condition, values_from = all_of(value_col))
    if (!all(c("PA", "AC") %in% colnames(wide))) return(NULL)
    wide <- wide[complete.cases(wide[, c("PA", "AC")]), , drop = FALSE]
    if (nrow(wide) == 0) return(NULL)
    diff <- wide$AC - wide$PA
    p <- if (length(diff) >= 2) {
      tryCatch(
        wilcox.test(wide$AC, wide$PA, paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    } else NA_real_
    dz <- if (length(diff) >= 2 && is.finite(sd(diff)) && sd(diff) > 0) {
      mean(diff) / sd(diff)
    } else NA_real_
    expected <- if (!is.null(expected_map) && feat %in% names(expected_map)) {
      expected_map[[feat]]
    } else NA_character_
    n_expected <- if (is.na(expected)) {
      NA_integer_
    } else if (expected == "AC>PA") {
      sum(diff > 0)
    } else {
      sum(diff < 0)
    }
    get_diff <- function(pid) {
      idx <- match(pid, wide$patient_id)
      if (is.na(idx)) NA_real_ else diff[idx]
    }
    data.frame(
      feature = feat,
      n_pairs = length(diff),
      mean_AC_minus_PA = mean(diff),
      median_AC_minus_PA = median(diff),
      paired_Cohen_dz = dz,
      wilcoxon_p = p,
      expected_direction = expected,
      n_pairs_in_expected_direction = n_expected,
      fraction_in_expected_direction = ifelse(is.na(n_expected), NA_real_, n_expected / length(diff)),
      P1_AC_minus_PA = get_diff("P1"),
      P2_AC_minus_PA = get_diff("P2"),
      P3_AC_minus_PA = get_diff("P3")
    )
  })
  bind_rows(out)
}

# =============================================================================
# A. 读取聚合矩阵并按后缀拆成 6 个样本
# =============================================================================
meta <- read.csv(file.path(BASE, "sample_metadata.csv"), check.names = FALSE)
meta$condition <- factor(meta$condition, levels = c("PA", "AC"))
write.csv(meta, file.path(OUT_DIR, "00_sample_metadata.csv"), row.names = FALSE)

agg_mat <- read_10x_any_compression(AGG_DIR)
barcode_suffix <- sub("^.*-([0-9]+)$", "\\1", colnames(agg_mat))
if (any(!grepl("^[1-6]$", barcode_suffix))) {
  stop("Unexpected GEM suffixes; expected only -1 through -6.")
}

suffix_counts <- as.data.frame(table(gem_suffix = barcode_suffix), stringsAsFactors = FALSE)
suffix_counts$gem_suffix <- as.integer(as.character(suffix_counts$gem_suffix))
colnames(suffix_counts)[2] <- "cells_in_processed_matrix"
suffix_counts <- left_join(meta, suffix_counts, by = "gem_suffix")
write.csv(suffix_counts, file.path(OUT_DIR, "00_barcode_suffix_sample_counts.csv"), row.names = FALSE)

obj_list <- list()
qc_rows <- list()
for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  suf <- as.character(meta$gem_suffix[i])
  cells <- colnames(agg_mat)[barcode_suffix == suf]
  if (length(cells) == 0) stop("No cells for suffix ", suf)

  mat <- agg_mat[, cells, drop = FALSE]
  colnames(mat) <- paste0(sid, "_", colnames(mat))
  obj <- CreateSeuratObject(
    counts = mat, project = "GSE159677",
    min.cells = 3, min.features = 200
  )
  obj$sample_id <- sid
  obj$patient_id <- meta$patient_id[i]
  obj$condition <- as.character(meta$condition[i])
  obj$tissue <- meta$tissue[i]
  obj$GSM <- meta$GSM[i]
  obj$gem_suffix <- meta$gem_suffix[i]
  obj$dataset_id <- "GSE159677"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  before_n <- ncol(obj)
  med_feat <- median(obj$nFeature_RNA)
  mad_feat <- mad(obj$nFeature_RNA, constant = 1.4826)
  med_count <- median(obj$nCount_RNA)
  mad_count <- mad(obj$nCount_RNA, constant = 1.4826)

  upper_feat <- max(3000, min(8000, med_feat + 4 * mad_feat))
  upper_count <- max(15000, min(60000, med_count + 4 * mad_count))
  keep <- obj$nFeature_RNA >= 200 &
    obj$nFeature_RNA <= upper_feat &
    obj$nCount_RNA <= upper_count &
    obj$percent.mt <= 20

  obj <- subset(obj, cells = colnames(obj)[keep])
  after_n <- ncol(obj)
  qc_rows[[sid]] <- data.frame(
    gem_suffix = meta$gem_suffix[i],
    sample_id = sid,
    patient_id = meta$patient_id[i],
    condition = as.character(meta$condition[i]),
    cells_before_qc = before_n,
    cells_after_qc = after_n,
    retain_rate = after_n / before_n,
    upper_nFeature = upper_feat,
    upper_nCount = upper_count,
    max_percent_mt = 20
  )
  obj_list[[sid]] <- obj
  rm(mat, obj)
  gc()
}
rm(agg_mat)
gc()

qc_summary <- bind_rows(qc_rows)
write.csv(qc_summary, file.path(OUT_DIR, "01_QC_summary_by_sample.csv"), row.names = FALSE)

qc_meta <- bind_rows(lapply(obj_list, function(x) {
  data.frame(
    sample_id = x$sample_id,
    condition = x$condition,
    nFeature_RNA = x$nFeature_RNA,
    nCount_RNA = x$nCount_RNA,
    percent.mt = x$percent.mt
  )
}))

p_qc1 <- ggplot(qc_meta, aes(sample_id, nFeature_RNA, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  labs(x = NULL, y = "Detected genes", title = "GSE159677 QC: genes per cell") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_qc2 <- ggplot(qc_meta, aes(sample_id, percent.mt, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  labs(x = NULL, y = "Mitochondrial percentage", title = "GSE159677 QC: mitochondrial fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png(p_qc1 / p_qc2, "01_QC_violin_plots.png", 12, 10)

# =============================================================================
# B. 全细胞图谱：PCA -> Harmony(patient_id) -> UMAP/聚类
# =============================================================================
full <- merge(x = obj_list[[1]], y = obj_list[-1], project = "GSE159677")
full <- join_layers_safe(full, "RNA")
DefaultAssay(full) <- "RNA"
full <- NormalizeData(full, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
full <- FindVariableFeatures(full, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
full <- ScaleData(full, features = VariableFeatures(full), verbose = FALSE)
full <- RunPCA(full, features = VariableFeatures(full), npcs = 40, verbose = FALSE)
full <- run_harmony_strict(full, group_var = "patient_id", reduction_name = "pca", save_name = "harmony")

harmony_dims_full <- seq_len(min(30, ncol(Embeddings(full, "harmony"))))
full <- RunUMAP(full, reduction = "harmony", dims = harmony_dims_full, seed.use = 20260731, verbose = FALSE)
full <- FindNeighbors(full, reduction = "harmony", dims = harmony_dims_full, verbose = FALSE)
full <- FindClusters(full, resolution = 0.5, random.seed = 20260731, verbose = FALSE)
full$full_cluster <- as.character(full$seurat_clusters)

p_full_condition <- DimPlot(full, reduction = "umap", group.by = "condition", pt.size = 0.15) +
  ggtitle("GSE159677 full atlas: AC versus PA")
p_full_sample <- DimPlot(full, reduction = "umap", group.by = "sample_id", pt.size = 0.12) +
  ggtitle("GSE159677 full atlas: Harmony corrected by patient")
p_full_cluster <- DimPlot(full, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.12) +
  ggtitle("GSE159677 full atlas: unsupervised clusters")
save_png(p_full_condition, "02_full_UMAP_condition_Harmony.png", 9, 7)
save_png(p_full_sample, "03_full_UMAP_sample_Harmony.png", 10, 7)
save_png(p_full_cluster, "04_full_UMAP_clusters_Harmony.png", 10, 8)

# =============================================================================
# C. cluster 级 + cell 级双重严格 EC 过滤
# =============================================================================
lineage_sets <- list(
  EC_core = c("PECAM1", "VWF", "CDH5", "KDR", "ESAM", "EMCN", "ENG", "RAMP2", "PLVAP", "EGFL7", "ERG", "EPAS1", "CLDN5", "TEK"),
  Immune = c("PTPRC", "LST1", "TYROBP", "FCER1G", "LYZ", "CD3D", "CD3E", "MS4A1", "CD79A", "NKG7"),
  SMC_Pericyte = c("ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "CSPG4", "MCAM", "PDGFRB", "DES"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PDGFRA", "C7", "COL6A1", "PI16", "CCDC80")
)

full_counts <- get_assay_matrix(full, "RNA", "counts")
full_data <- get_assay_matrix(full, "RNA", "data")

full$n_EC_markers <- marker_count(full_counts, lineage_sets$EC_core)
full$n_immune_markers <- marker_count(full_counts, lineage_sets$Immune)
full$n_SMC_markers <- marker_count(full_counts, lineage_sets$SMC_Pericyte)
full$n_fibro_markers <- marker_count(full_counts, lineage_sets$Fibroblast)
full$PTPRC_detected <- if ("PTPRC" %in% rownames(full_counts)) {
  as.numeric(full_counts["PTPRC", ] > 0)
} else 0

full$EC_core_z <- gene_z_score(full_data, lineage_sets$EC_core)
full$Immune_z <- gene_z_score(full_data, lineage_sets$Immune)
full$SMC_Pericyte_z <- gene_z_score(full_data, lineage_sets$SMC_Pericyte)
full$Fibroblast_z <- gene_z_score(full_data, lineage_sets$Fibroblast)

cluster_summary <- full@meta.data %>%
  group_by(full_cluster) %>%
  summarise(
    n_cells = n(),
    frac_EC_2plus = mean(n_EC_markers >= 2),
    frac_EC_3plus = mean(n_EC_markers >= 3),
    frac_PTPRC_positive = mean(PTPRC_detected > 0),
    frac_strong_immune = mean(n_immune_markers >= 2),
    frac_strong_SMC = mean(n_SMC_markers >= 2),
    frac_strong_fibro = mean(n_fibro_markers >= 2),
    median_EC_z = median(EC_core_z, na.rm = TRUE),
    median_Immune_z = median(Immune_z, na.rm = TRUE),
    median_SMC_z = median(SMC_Pericyte_z, na.rm = TRUE),
    median_Fibro_z = median(Fibroblast_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median_max_nonEC_z = pmax(median_Immune_z, median_SMC_z, median_Fibro_z),
    strict_EC_cluster =
      frac_EC_2plus >= 0.35 &
      median_EC_z >= median_max_nonEC_z + 0.10 &
      frac_PTPRC_positive <= 0.10 &
      frac_strong_immune <= 0.15 &
      frac_strong_SMC <= 0.20 &
      frac_strong_fibro <= 0.20,
    broad_EC_cluster =
      frac_EC_2plus >= 0.20 &
      median_EC_z > median_Immune_z
  )
write.csv(cluster_summary, file.path(OUT_DIR, "02_cluster_purity_audit_STRICT.csv"), row.names = FALSE)

strict_clusters <- cluster_summary$full_cluster[cluster_summary$strict_EC_cluster]
if (length(strict_clusters) == 0) {
  stop("No strict EC clusters passed. Review 02_cluster_purity_audit_STRICT.csv")
}
message("Strict EC clusters: ", paste(strict_clusters, collapse = ", "))

full$strict_EC_cell <-
  full$full_cluster %in% strict_clusters &
  full$n_EC_markers >= 2 &
  full$PTPRC_detected == 0 &
  full$n_immune_markers <= 1 &
  full$n_SMC_markers <= 1 &
  full$n_fibro_markers <= 1 &
  full$EC_core_z >= pmax(full$Immune_z, full$SMC_Pericyte_z, full$Fibroblast_z) + 0.05

strict_cells <- rownames(full@meta.data)[full$strict_EC_cell]
if (length(strict_cells) < 500) {
  stop("Strict filter retained only ", length(strict_cells), " cells; manual review is required.")
}

strict_count_table <- full@meta.data %>%
  count(sample_id, patient_id, condition, strict_EC_cell, name = "n_cells") %>%
  arrange(sample_id, desc(strict_EC_cell))
write.csv(strict_count_table, file.path(OUT_DIR, "03_strict_EC_filter_counts.csv"), row.names = FALSE)

strict_by_sample <- full@meta.data[full$strict_EC_cell, , drop = FALSE] %>%
  count(sample_id, patient_id, condition, name = "n_strict_EC")
write.csv(strict_by_sample, file.path(OUT_DIR, "04_strict_EC_counts_by_sample.csv"), row.names = FALSE)

missing_samples <- setdiff(meta$sample_id, strict_by_sample$sample_id)
if (length(missing_samples) > 0) stop("No strict EC retained in: ", paste(missing_samples, collapse = ", "))
if (any(strict_by_sample$n_strict_EC < 30)) {
  warning("At least one sample has fewer than 30 strict EC; inspect 04_strict_EC_counts_by_sample.csv")
}

marker_panel_full <- unique(c(
  lineage_sets$EC_core,
  lineage_sets$Immune,
  lineage_sets$SMC_Pericyte,
  lineage_sets$Fibroblast
))
marker_panel_full <- intersect(marker_panel_full, rownames(full))
p_dot_full <- DotPlot(full, features = marker_panel_full, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("Full-atlas lineage audit before strict EC filtering")
save_png(p_dot_full, "05_full_cluster_marker_dotplot_STRICT.png", 18, 9)

p_filter <- DimPlot(full, group.by = "strict_EC_cell", pt.size = 0.12) +
  ggtitle("High-confidence EC retained by strict cluster + cell filtering")
save_png(p_filter, "06_full_UMAP_strict_EC_selection.png", 10, 8)

# =============================================================================
# D. 严格 EC 重聚类：Harmony 失败即停止
# =============================================================================
ec <- subset(full, cells = strict_cells)
ec <- join_layers_safe(ec, "RNA")
DefaultAssay(ec) <- "RNA"
ec <- NormalizeData(ec, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ec <- FindVariableFeatures(ec, selection.method = "vst", nfeatures = min(2500, nrow(ec)), verbose = FALSE)
ec <- ScaleData(ec, features = VariableFeatures(ec), verbose = FALSE)
ec <- RunPCA(ec, features = VariableFeatures(ec), npcs = 40, verbose = FALSE)
ec <- run_harmony_strict(ec, group_var = "patient_id", reduction_name = "pca", save_name = "harmony")

harmony_dims_ec <- seq_len(min(30, ncol(Embeddings(ec, "harmony"))))
ec <- RunUMAP(ec, reduction = "harmony", dims = harmony_dims_ec, seed.use = 20260731, verbose = FALSE)
ec <- FindNeighbors(ec, reduction = "harmony", dims = harmony_dims_ec, verbose = FALSE)
ec <- FindClusters(ec, resolution = 0.4, random.seed = 20260731, verbose = FALSE)
ec$strict_EC_cluster <- as.character(ec$seurat_clusters)

p_ec_condition <- DimPlot(ec, group.by = "condition", pt.size = 0.30) +
  ggtitle("Strict EC: AC versus PA")
p_ec_sample <- DimPlot(ec, group.by = "sample_id", pt.size = 0.25) +
  ggtitle("Strict EC: sample distribution after patient-level Harmony")
p_ec_cluster <- DimPlot(ec, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.25) +
  ggtitle("Strict EC reclustering")
save_png(p_ec_condition, "07_strict_EC_UMAP_condition.png", 9, 7)
save_png(p_ec_sample, "08_strict_EC_UMAP_sample_Harmony.png", 10, 7)
save_png(p_ec_cluster, "09_strict_EC_UMAP_clusters.png", 9, 7)

# 严格 EC 纯度复核
strict_marker_panel <- intersect(unique(c(
  "PECAM1", "VWF", "CDH5", "KDR", "ESAM", "CLDN5", "RAMP2", "ERG",
  "PTPRC", "LST1", "TYROBP",
  "ACTA2", "TAGLN", "MYH11", "RGS5", "PDGFRB",
  "COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PI16"
)), rownames(ec))
p_strict_dot <- DotPlot(ec, features = strict_marker_panel, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("Purity audit after strict EC filtering")
save_png(p_strict_dot, "10_strict_EC_purity_marker_dotplot.png", 16, 8)

# =============================================================================
# E. 固定验证程序：不再把动脉身份等同于稳态
# =============================================================================
signatures <- list(
  C7_arterial_identity = c("EFNB2", "GJA4", "GJA5", "SOX17", "HEY1", "SEMA3G", "CXCL12"),
  C7_homeostatic_protection = c("KLF2", "KLF4", "NOS3", "THBD", "TEK", "ERG", "CLDN5", "RAMP2", "EPAS1"),
  C3_ANGPT2_PGF_activation = c("ANGPT2", "PGF"),
  C1_C5_inflammatory_adhesion = c("VCAM1", "SELE", "ICAM1", "CXCL2", "CXCL3", "NFKBIA", "IFITM3", "IFIT1", "ISG15", "STAT1"),
  C2_vascular_remodeling = c("ITLN1", "SULF1", "MGP", "CLU", "ELN", "COL8A1", "COL8A2", "EFEMP1", "MFAP4"),
  C4_ECM_remodeling_primary = c("POSTN", "IGFBP5", "SPARC", "COL6A2", "TIMP1", "MAFB", "NR2F2", "FN1", "VCAN"),
  C4_broad_ECM_sensitivity = c("POSTN", "IGFBP5", "SPARC", "COL6A2", "TIMP1", "MAFB", "NR2F2", "FN1", "VCAN", "COL1A1", "COL1A2", "COL3A1"),
  Fibroblast_contamination_QC = lineage_sets$Fibroblast,
  SMC_Pericyte_contamination_QC = lineage_sets$SMC_Pericyte
)

signature_role <- c(
  C7_arterial_identity = "identity_no_fixed_disease_direction",
  C7_homeostatic_protection = "primary_validation",
  C3_ANGPT2_PGF_activation = "primary_validation",
  C1_C5_inflammatory_adhesion = "primary_validation",
  C2_vascular_remodeling = "primary_validation",
  C4_ECM_remodeling_primary = "primary_validation",
  C4_broad_ECM_sensitivity = "sensitivity_analysis",
  Fibroblast_contamination_QC = "purity_QC",
  SMC_Pericyte_contamination_QC = "purity_QC"
)

expected_program <- c(
  C7_arterial_identity = NA_character_,
  C7_homeostatic_protection = "PA>AC",
  C3_ANGPT2_PGF_activation = "AC>PA",
  C1_C5_inflammatory_adhesion = "AC>PA",
  C2_vascular_remodeling = "AC>PA",
  C4_ECM_remodeling_primary = "AC>PA",
  C4_broad_ECM_sensitivity = "AC>PA",
  Fibroblast_contamination_QC = NA_character_,
  SMC_Pericyte_contamination_QC = NA_character_
)

signatures_available <- lapply(signatures, function(x) unique(intersect(x, rownames(ec))))
sig_table <- bind_rows(lapply(names(signatures_available), function(nm) {
  data.frame(
    program = nm,
    role = signature_role[[nm]],
    expected_direction = expected_program[[nm]],
    gene = signatures_available[[nm]]
  )
}))
write.csv(sig_table, file.path(OUT_DIR, "05_fixed_validation_signatures_STRICT.csv"), row.names = FALSE)

short_sig <- names(signatures_available)[lengths(signatures_available) < 3]
if (length(short_sig) > 0) {
  warning("Signatures with fewer than 3 detected genes: ", paste(short_sig, collapse = ", "))
}

ec_data <- get_assay_matrix(ec, "RNA", "data")
for (nm in names(signatures_available)) {
  ec[[paste0(nm, "_meanlog")]] <- mean_log_score(ec_data, signatures_available[[nm]])
  ec[[paste0(nm, "_z")]] <- gene_z_score(ec_data, signatures_available[[nm]])
}

# =============================================================================
# F. 参考对象标签转移；低置信度标为 Uncertain
# =============================================================================
normalize_cluster_label <- function(x) {
  x <- as.character(x)
  x <- gsub("^cluster[_ ]*", "", x, ignore.case = TRUE)
  x <- gsub("^C", "", x, ignore.case = TRUE)
  x
}

find_reference_cluster_labels <- function(ref) {
  wanted <- c("1", "2", "3", "4", "5", "7")

  valid_cluster_vector <- function(v) {
    vv <- normalize_cluster_label(v)
    uu <- unique(vv[!is.na(vv) & vv != ""])
    length(uu) >= length(wanted) &&
      length(uu) <= 40 &&
      all(grepl("^[0-9]+$", uu)) &&
      all(wanted %in% uu)
  }

  idv <- normalize_cluster_label(Idents(ref))
  if (valid_cluster_vector(idv)) return(idv)

  preferred <- c(
    "seurat_clusters", "cluster", "clusters", "cluster_id",
    "EC_cluster", "ec_cluster", "sling_cluster", "slingshot_cluster"
  )
  candidates <- unique(c(
    preferred[preferred %in% colnames(ref@meta.data)],
    colnames(ref@meta.data)
  ))

  for (cn in candidates) {
    vv <- ref@meta.data[[cn]]
    if (valid_cluster_vector(vv)) {
      message("Using reference cluster metadata column: ", cn)
      return(normalize_cluster_label(vv))
    }
  }
  stop("Could not safely identify reference labels 1,2,3,4,5,7")
}

if (!file.exists(REF_EC)) {
  hits <- list.files(
    "/root", pattern = "slingshot.*scored\\.rds$|arterial_dysfunction.*\\.rds$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  if (length(hits) == 0) stop("Discovery EC reference RDS was not found.")
  REF_EC <- hits[1]
}
message("Loading discovery EC reference: ", REF_EC)
ref <- readRDS(REF_EC)
if (!"RNA" %in% Assays(ref)) stop("Reference object has no RNA assay.")
ref <- join_layers_safe(ref, "RNA")
DefaultAssay(ref) <- "RNA"
ref$reference_cluster_num <- find_reference_cluster_labels(ref)
ref <- subset(ref, subset = reference_cluster_num %in% c("1", "2", "3", "4", "5", "7"))
ref$reference_state <- paste0("C", ref$reference_cluster_num)
Idents(ref) <- ref$reference_state
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = min(3000, nrow(ref)), verbose = FALSE)
ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 30, verbose = FALSE)

common_features <- intersect(VariableFeatures(ref), rownames(ec))
if (length(common_features) < 500) {
  stop("Only ", length(common_features), " common variable genes; reference mapping is unreliable.")
}

anchors <- FindTransferAnchors(
  reference = ref,
  query = ec,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  features = common_features,
  dims = 1:30,
  verbose = FALSE
)
pred <- TransferData(
  anchorset = anchors,
  refdata = ref$reference_state,
  dims = 1:30,
  verbose = FALSE
)

score_cols <- grep("^prediction\\.score\\.C", colnames(pred), value = TRUE)
if (length(score_cols) >= 2) {
  score_matrix <- as.matrix(pred[, score_cols, drop = FALSE])
  prediction_margin <- apply(score_matrix, 1, function(x) {
    sx <- sort(x, decreasing = TRUE)
    sx[1] - sx[2]
  })
} else {
  prediction_margin <- rep(NA_real_, nrow(pred))
}
pred$prediction.margin <- prediction_margin
ec <- AddMetaData(ec, pred)
ec$mapped_state <- ifelse(
  ec$prediction.score.max >= 0.50 &
    (is.na(ec$prediction.margin) | ec$prediction.margin >= 0.05),
  as.character(ec$predicted.id),
  "Uncertain"
)

all_states <- c("C1", "C2", "C3", "C4", "C5", "C7", "Uncertain")
map_counts <- ec@meta.data %>%
  count(sample_id, mapped_state, name = "n_cells")
sample_totals <- ec@meta.data %>%
  count(sample_id, name = "n_total_strict_EC")
map_diag <- tidyr::expand_grid(
  sample_id = meta$sample_id,
  mapped_state = all_states
) %>%
  left_join(map_counts, by = c("sample_id", "mapped_state")) %>%
  mutate(n_cells = tidyr::replace_na(n_cells, 0L)) %>%
  left_join(
    meta %>% select(sample_id, patient_id, condition),
    by = "sample_id"
  ) %>%
  left_join(sample_totals, by = "sample_id") %>%
  mutate(proportion = n_cells / n_total_strict_EC)
write.csv(map_diag, file.path(OUT_DIR, "06_mapped_state_counts_proportions_STRICT.csv"), row.names = FALSE)

p_state <- DimPlot(ec, group.by = "mapped_state", pt.size = 0.28) +
  ggtitle("Strict EC mapped to discovery states")
p_conf <- FeaturePlot(ec, features = "prediction.score.max", pt.size = 0.28) +
  ggtitle("Reference mapping confidence")
save_png(p_state, "11_strict_EC_UMAP_mapped_state.png", 10, 7)
save_png(p_conf, "12_strict_EC_mapping_confidence.png", 9, 7)

state_marker_panel <- intersect(unique(c(
  "PECAM1", "VWF", "CDH5", "KDR", "RAMP2",
  signatures$C7_arterial_identity,
  signatures$C7_homeostatic_protection,
  signatures$C3_ANGPT2_PGF_activation,
  "CCL14",
  signatures$C1_C5_inflammatory_adhesion,
  signatures$C2_vascular_remodeling,
  signatures$C4_ECM_remodeling_primary,
  "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"
)), rownames(ec))
p_state_dot <- DotPlot(ec, features = state_marker_panel, group.by = "mapped_state") +
  RotatedAxis() + ggtitle("Fixed discovery markers in strict GSE159677 EC")
save_png(p_state_dot, "13_strict_EC_state_marker_dotplot.png", 20, 9)

# =============================================================================
# G. 样本级固定程序配对验证
# =============================================================================
program_names <- names(signatures_available)
program_z_cols <- paste0(program_names, "_z")

sample_program <- ec@meta.data %>%
  group_by(sample_id, patient_id, condition) %>%
  summarise(
    n_strict_EC = n(),
    across(all_of(program_z_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
colnames(sample_program)[match(program_z_cols, colnames(sample_program))] <- program_names
write.csv(sample_program, file.path(OUT_DIR, "07_sample_program_scores_STRICT.csv"), row.names = FALSE)

program_long <- sample_program %>%
  pivot_longer(cols = all_of(program_names), names_to = "program", values_to = "mean_z_score")

p_program <- ggplot(
  program_long,
  aes(condition, mean_z_score, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 2.6) +
  facet_wrap(~ program, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Mean standardized program score", title = "Strict EC paired validation: PA versus AC")
save_png(p_program, "14_paired_program_scores_STRICT.png", 13, 14)

primary_heat_programs <- names(signature_role)[signature_role %in% c(
  "primary_validation", "identity_no_fixed_disease_direction", "sensitivity_analysis"
)]
heat_mat <- as.data.frame(sample_program)
rownames(heat_mat) <- heat_mat$sample_id
heat_mat <- as.matrix(heat_mat[, primary_heat_programs, drop = FALSE])
heat_z <- t(scale(t(heat_mat)))
heat_z[!is.finite(heat_z)] <- 0
png(file.path(OUT_DIR, "15_sample_program_heatmap_STRICT.png"), width = 3200, height = 2300, res = 300, bg = "white")
pheatmap::pheatmap(
  heat_z, cluster_rows = TRUE, cluster_cols = FALSE,
  main = "Strict EC program scores (row z-score)", border_color = NA
)
dev.off()

program_stats <- paired_stats(program_long, "mean_z_score", "program", expected_program)
program_stats <- left_join(
  program_stats,
  data.frame(program = names(signature_role), role = unname(signature_role)),
  by = c("feature" = "program")
)
write.csv(program_stats, file.path(OUT_DIR, "08_paired_program_statistics_STRICT.csv"), row.names = FALSE)

# =============================================================================
# H. 映射状态比例的配对比较
# =============================================================================
p_state_pair <- ggplot(
  map_diag,
  aes(condition, proportion, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 2.5) +
  facet_wrap(~ mapped_state, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Proportion among strict EC", title = "Strict EC mapped-state proportions")
save_png(p_state_pair, "16_paired_mapped_state_proportions_STRICT.png", 13, 9)

state_expected <- c(
  C7 = "PA>AC", C3 = "AC>PA", C1 = "AC>PA", C5 = "AC>PA",
  C2 = "AC>PA", C4 = "AC>PA", Uncertain = NA_character_
)
state_stats <- paired_stats(map_diag, "proportion", "mapped_state", state_expected)
write.csv(state_stats, file.path(OUT_DIR, "09_paired_state_proportion_statistics_STRICT.csv"), row.names = FALSE)

# =============================================================================
# I. 关键基因配对表达：保留原假设，同时把 CCL14 单独报告
# =============================================================================
key_gene_groups <- list(
  arterial_identity = signatures$C7_arterial_identity,
  homeostatic_protection = signatures$C7_homeostatic_protection,
  C3_activation = unique(c(signatures$C3_ANGPT2_PGF_activation, "CCL14")),
  inflammatory_adhesion = signatures$C1_C5_inflammatory_adhesion,
  C2_remodeling = signatures$C2_vascular_remodeling,
  C4_ECM = signatures$C4_ECM_remodeling_primary
)
key_genes <- intersect(unique(unlist(key_gene_groups)), rownames(ec))

key_rows <- lapply(unique(ec$sample_id), function(sid) {
  cells <- colnames(ec)[ec$sample_id == sid]
  md1 <- ec@meta.data[cells[1], c("sample_id", "patient_id", "condition"), drop = FALSE]
  vals <- Matrix::rowMeans(ec_data[key_genes, cells, drop = FALSE])
  data.frame(
    sample_id = md1$sample_id,
    patient_id = md1$patient_id,
    condition = md1$condition,
    gene = names(vals),
    mean_log_normalized_expression = as.numeric(vals)
  )
})
key_gene_sample <- bind_rows(key_rows)
write.csv(key_gene_sample, file.path(OUT_DIR, "10_key_gene_sample_expression_STRICT.csv"), row.names = FALSE)

p_key <- ggplot(
  key_gene_sample,
  aes(condition, mean_log_normalized_expression, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = "Mean log-normalized expression", title = "Strict EC key genes: paired PA versus AC")
save_png(p_key, "17_key_gene_paired_expression_STRICT.png", 17, 15)

key_expected <- setNames(rep(NA_character_, length(key_genes)), key_genes)
key_expected[intersect(names(key_expected), signatures$C7_homeostatic_protection)] <- "PA>AC"
key_expected[intersect(names(key_expected), unique(c(
  signatures$C3_ANGPT2_PGF_activation, "CCL14",
  signatures$C1_C5_inflammatory_adhesion,
  signatures$C2_vascular_remodeling,
  signatures$C4_ECM_remodeling_primary
)))] <- "AC>PA"
key_stats <- paired_stats(
  key_gene_sample,
  "mean_log_normalized_expression",
  "gene",
  key_expected
)
write.csv(key_stats, file.path(OUT_DIR, "11_key_gene_paired_statistics_STRICT.csv"), row.names = FALSE)

# =============================================================================
# J. C4 污染鲁棒性：进一步去掉任何 fibro/SMC marker 阳性的 ultra-strict EC
# =============================================================================
ultra_cells <- rownames(ec@meta.data)[
  ec$n_EC_markers >= 3 &
    ec$n_immune_markers == 0 &
    ec$n_SMC_markers == 0 &
    ec$n_fibro_markers == 0
]
ultra_count <- ec@meta.data[ultra_cells, , drop = FALSE] %>%
  count(sample_id, patient_id, condition, name = "n_ultra_strict_EC")
write.csv(ultra_count, file.path(OUT_DIR, "12_ultra_strict_EC_counts.csv"), row.names = FALSE)

c4_ultra_stats <- data.frame()
if (length(ultra_cells) >= 200 && all(meta$sample_id %in% ultra_count$sample_id)) {
  c4_cols <- paste0(c("C4_ECM_remodeling_primary", "C4_broad_ECM_sensitivity"), "_z")
  c4_ultra <- ec@meta.data[ultra_cells, , drop = FALSE] %>%
    group_by(sample_id, patient_id, condition) %>%
    summarise(across(all_of(c4_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  colnames(c4_ultra)[match(c4_cols, colnames(c4_ultra))] <- c(
    "C4_ECM_remodeling_primary", "C4_broad_ECM_sensitivity"
  )
  write.csv(c4_ultra, file.path(OUT_DIR, "13_C4_ultra_strict_sample_scores.csv"), row.names = FALSE)

  c4_ultra_long <- c4_ultra %>%
    pivot_longer(
      cols = c("C4_ECM_remodeling_primary", "C4_broad_ECM_sensitivity"),
      names_to = "program", values_to = "mean_z_score"
    )
  c4_ultra_stats <- paired_stats(c4_ultra_long, "mean_z_score", "program", expected_program)
  write.csv(c4_ultra_stats, file.path(OUT_DIR, "14_C4_ultra_strict_paired_statistics.csv"), row.names = FALSE)

  p_c4_ultra <- ggplot(
    c4_ultra_long,
    aes(condition, mean_z_score, group = patient_id, colour = patient_id)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~ program, scales = "free_y") +
    labs(x = NULL, y = "Mean standardized score", title = "C4 robustness after ultra-strict EC filtering")
  save_png(p_c4_ultra, "18_C4_ultra_strict_robustness.png", 10, 6)
} else {
  writeLines(
    "Ultra-strict C4 sensitivity analysis was skipped because too few cells or incomplete sample coverage remained.",
    file.path(OUT_DIR, "14_C4_ultra_strict_NOT_RUN.txt")
  )
}

c4_cor <- data.frame(
  comparison = c(
    "C4_primary_vs_fibroblast_QC",
    "C4_broad_vs_fibroblast_QC",
    "C4_primary_vs_SMC_QC"
  ),
  spearman_rho = c(
    suppressWarnings(cor(ec$C4_ECM_remodeling_primary_z, ec$Fibroblast_contamination_QC_z, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(ec$C4_broad_ECM_sensitivity_z, ec$Fibroblast_contamination_QC_z, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(ec$C4_ECM_remodeling_primary_z, ec$SMC_Pericyte_contamination_QC_z, method = "spearman", use = "complete.obs"))
  )
)
write.csv(c4_cor, file.path(OUT_DIR, "15_C4_contamination_score_correlations.csv"), row.names = FALSE)

# =============================================================================
# K. strict EC pseudobulk paired edgeR
# =============================================================================
ec_counts <- get_assay_matrix(ec, "RNA", "counts")
samples <- meta$sample_id
pb_counts <- sapply(samples, function(sid) {
  cells <- colnames(ec)[ec$sample_id == sid]
  if (length(cells) == 0) stop("No strict EC for pseudobulk sample ", sid)
  Matrix::rowSums(ec_counts[, cells, drop = FALSE])
})
if (!is.matrix(pb_counts)) pb_counts <- as.matrix(pb_counts)
colnames(pb_counts) <- samples
write.csv(
  data.frame(gene = rownames(pb_counts), pb_counts, check.names = FALSE),
  file.path(OBJ_DIR, "strict_EC_pseudobulk_raw_counts.csv"),
  row.names = FALSE
)

if (!requireNamespace("edgeR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    try(install.packages("BiocManager", repos = "https://cloud.r-project.org"), silent = TRUE)
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    try(BiocManager::install("edgeR", ask = FALSE, update = FALSE), silent = TRUE)
  }
}

edgeR_success <- FALSE
if (requireNamespace("edgeR", quietly = TRUE)) {
  edgeR_success <- tryCatch({
    suppressPackageStartupMessages(library(edgeR))
    pb_meta <- meta[match(colnames(pb_counts), meta$sample_id), ]
    pb_meta$condition <- factor(pb_meta$condition, levels = c("PA", "AC"))
    pb_meta$patient_id <- factor(pb_meta$patient_id)
    design <- model.matrix(~ patient_id + condition, data = pb_meta)
    if (qr(design)$rank < ncol(design)) stop("Paired design matrix is not full rank.")

    y <- edgeR::DGEList(counts = pb_counts)
    keep_gene <- edgeR::filterByExpr(y, design = design, min.count = 5)
    y <- y[keep_gene, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    coef_name <- grep("conditionAC", colnames(design), value = TRUE)
    if (length(coef_name) != 1) stop("conditionAC coefficient not found.")
    qlf <- edgeR::glmQLFTest(fit, coef = which(colnames(design) == coef_name))
    tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
    tab$gene <- rownames(tab)
    tab <- tab[, c("gene", setdiff(colnames(tab), "gene"))]
    write.csv(tab, file.path(OBJ_DIR, "strict_EC_edgeR_all_genes.csv"), row.names = FALSE)
    write.csv(
      tab[tab$gene %in% key_genes, , drop = FALSE],
      file.path(OUT_DIR, "16_strict_EC_edgeR_key_genes.csv"),
      row.names = FALSE
    )
    TRUE
  }, error = function(e) {
    writeLines(
      paste("edgeR paired analysis failed:", conditionMessage(e)),
      file.path(OUT_DIR, "16_edgeR_NOT_RUN.txt")
    )
    FALSE
  })
} else {
  writeLines(
    "edgeR was unavailable. Raw strict-EC pseudobulk counts were saved, but formal testing was skipped.",
    file.path(OUT_DIR, "16_edgeR_NOT_RUN.txt")
  )
}

# =============================================================================
# L. 汇总结论等级、保存对象和会话信息
# =============================================================================
primary_programs <- c(
  "C7_homeostatic_protection",
  "C3_ANGPT2_PGF_activation",
  "C1_C5_inflammatory_adhesion",
  "C2_vascular_remodeling",
  "C4_ECM_remodeling_primary"
)
summary_program <- program_stats %>%
  filter(feature %in% primary_programs) %>%
  mutate(
    evidence_grade = case_when(
      n_pairs_in_expected_direction == 3 ~ "directionally_consistent_3_of_3",
      n_pairs_in_expected_direction == 2 ~ "partially_consistent_2_of_3",
      n_pairs_in_expected_direction <= 1 ~ "not_consistent",
      TRUE ~ "not_assessed"
    )
  ) %>%
  select(
    validation_item = feature,
    expected_direction,
    mean_AC_minus_PA,
    paired_Cohen_dz,
    n_pairs_in_expected_direction,
    fraction_in_expected_direction,
    evidence_grade
  )
write.csv(summary_program, file.path(OUT_DIR, "17_primary_validation_summary_STRICT.csv"), row.names = FALSE)

saveRDS(full, file.path(OBJ_DIR, "GSE159677_full_Harmony_strictECflag.rds"), compress = "gzip")
saveRDS(ec, file.path(OBJ_DIR, "GSE159677_strict_EC_external_validation.rds"), compress = "gzip")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_STRICT.txt"))

summary_text <- c(
  "GSE159677 strict independent external validation completed.",
  paste0("Full cells after QC: ", ncol(full)),
  paste0("Strict EC cells: ", ncol(ec)),
  paste0("Strict full-atlas EC clusters: ", paste(strict_clusters, collapse = ", ")),
  paste0("Harmony full atlas by patient_id: ", "harmony" %in% Reductions(full)),
  paste0("Harmony strict EC by patient_id: ", "harmony" %in% Reductions(ec)),
  paste0("Reference mapping uncertain fraction: ", round(mean(ec$mapped_state == "Uncertain"), 4)),
  paste0("edgeR completed: ", edgeR_success),
  "Primary evidence is based on strict EC, patient-paired direction, effect size, and pseudobulk results.",
  "C7 arterial identity is reported separately and has no prespecified PA-versus-AC direction.",
  "C4 primary score excludes COL1A1/COL1A2/COL3A1; broad ECM and ultra-strict analyses are sensitivity checks.",
  "PA is proximal adjacent diseased artery, not a truly healthy control."
)
writeLines(summary_text, file.path(OUT_DIR, "README_STRICT_results.txt"))
cat(paste(summary_text, collapse = "\n"), "\n")
RSCRIPT

# -----------------------------------------------------------------------------
# 5. 运行 R
# -----------------------------------------------------------------------------
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: 当前环境没有 Rscript。请先进入之前运行 Seurat 的 R 环境。" >&2
  exit 1
fi

Rscript "$R_SCRIPT"

# -----------------------------------------------------------------------------
# 6. 自动生成轻量审阅包和大型对象包
# -----------------------------------------------------------------------------
REVIEW_DIR="${BASE}/GSE159677_strict_review_small"
rm -rf "$REVIEW_DIR"
mkdir -p "$REVIEW_DIR/results" "$REVIEW_DIR/logs"

find "$OUT_DIR" -maxdepth 1 -type f \
  \( -name '*.png' -o -name '*.csv' -o -name '*.txt' \) \
  -exec cp -f {} "$REVIEW_DIR/results/" \;

cp -f "${BASE}/sample_metadata.csv" "$REVIEW_DIR/"
cp -f "$R_SCRIPT" "$REVIEW_DIR/"
find "$LOG_DIR" -maxdepth 1 -type f -name '*.log' -exec cp -f {} "$REVIEW_DIR/logs/" \;
find "$REVIEW_DIR" -type f -printf '%P\t%k KB\n' | sort > "$REVIEW_DIR/file_manifest.txt"

tar -czf "${BASE}/GSE159677_strict_review_small.tar.gz" \
  -C "$BASE" "$(basename "$REVIEW_DIR")"

# 大型 RDS 和全基因表保留在 objects 目录，不自动再次压缩，避免耗时和超大压缩包。
cat <<EOF
============================================================
STRICT ANALYSIS COMPLETED: $(date)

请下载并上传给我的小型结果包：
  ${BASE}/GSE159677_strict_review_small.tar.gz

大型对象和全基因结果仅保留在服务器，暂时不要上传：
  ${BASE}/objects/

核心结果目录：
  ${OUT_DIR}

优先检查：
  ${OUT_DIR}/README_STRICT_results.txt
  ${OUT_DIR}/02_cluster_purity_audit_STRICT.csv
  ${OUT_DIR}/04_strict_EC_counts_by_sample.csv
  ${OUT_DIR}/08_paired_program_statistics_STRICT.csv
  ${OUT_DIR}/11_key_gene_paired_statistics_STRICT.csv
  ${OUT_DIR}/14_C4_ultra_strict_paired_statistics.csv
  ${OUT_DIR}/16_strict_EC_edgeR_key_genes.csv
  ${OUT_DIR}/17_primary_validation_summary_STRICT.csv
  ${OUT_DIR}/08_strict_EC_UMAP_sample_Harmony.png
  ${OUT_DIR}/10_strict_EC_purity_marker_dotplot.png
  ${OUT_DIR}/13_strict_EC_state_marker_dotplot.png
  ${OUT_DIR}/14_paired_program_scores_STRICT.png
  ${OUT_DIR}/18_C4_ultra_strict_robustness.png
============================================================
EOF
