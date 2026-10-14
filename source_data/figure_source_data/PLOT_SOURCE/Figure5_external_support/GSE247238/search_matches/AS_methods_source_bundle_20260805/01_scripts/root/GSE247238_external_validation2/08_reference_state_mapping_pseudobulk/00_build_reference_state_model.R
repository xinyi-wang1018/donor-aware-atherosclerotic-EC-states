suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
set.seed(20260801)

ref_path <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
outdir <- "/root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk/01_reference_model"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

TOP_MARKERS_PER_STATE <- 125
MIN_PCT <- 0.10
MIN_AVG_LOGEXPR <- 0.10

candidate_genes <- c(
  "TGFBR3", "SMAD1", "SOX4", "HIF1A", "ANGPT2",
  "NAMPT", "THBS1", "CTHRC1", "LOX", "SULF1"
)

# 候选和SULF1全部从分类特征中排除，避免循环验证。
excluded_classifier_genes <- unique(toupper(c(
  candidate_genes,
  grep("^(MT-|RPL|RPS|MIR|SNOR|HIST|HBA|HBB)", character(), value = TRUE)
)))

manual_anchors <- list(
  C7_quiescent = c(
    "CXCL12", "EFNB2", "GJA4", "GJA5", "SOX17", "HEY1",
    "KLF2", "KLF4", "THBD", "EPAS1", "BMX", "DLL4"
  ),
  C3_early_activation = c(
    "CCL14", "PGF", "KDR", "PLVAP", "VWF", "ESM1", "ENG"
  ),
  C5_transition = c(
    "ISG15", "IFI44L", "IFIT1", "IFIT3", "IFITM1", "IFITM3",
    "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CXCL10", "STAT1", "IRF1"
  ),
  C1C2_late_remodeling = c(
    "ITLN1", "MGP", "CLU", "ELN", "COL8A1", "EDN1", "PROCR",
    "FN1", "BGN", "FBLN5", "CCN1", "TIMP1"
  ),
  C4_independent_ECM = c(
    "POSTN", "IGFBP5", "SPARC", "COL6A2", "MAFB", "NR2F2",
    "COL15A1", "ACKR1", "PLVAP"
  )
)
manual_anchors <- lapply(
  manual_anchors,
  function(x) setdiff(unique(toupper(x)), excluded_classifier_genes)
)

standardize_cluster <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^[Cc]", "", x)
  ifelse(grepl("^[0-9]+$", x), paste0("C", x), x)
}

get_assay_matrix_safe <- function(object, layer_name) {
  tryCatch(
    GetAssayData(object, assay = "RNA", layer = layer_name),
    error = function(e1) {
      tryCatch(
        GetAssayData(object, assay = "RNA", slot = layer_name),
        error = function(e2) NULL
      )
    }
  )
}

message("读取参考对象...")
obj <- readRDS(ref_path)
if (!inherits(obj, "Seurat")) stop("参考输入不是Seurat对象。")
if (!("RNA" %in% names(obj@assays))) stop("参考对象没有RNA assay。")

cluster_candidates <- c("EC_cluster", "EC_cluster_std", "seurat_clusters")
cluster_col <- cluster_candidates[cluster_candidates %in% colnames(obj@meta.data)][1]
if (is.na(cluster_col)) stop("参考对象中未找到EC_cluster/EC_cluster_std/seurat_clusters列。")

obj$cluster_for_mapping <- standardize_cluster(obj@meta.data[[cluster_col]])
keep_clusters <- c("C7", "C3", "C5", "C1", "C2", "C4")
keep_cells <- rownames(obj@meta.data)[obj$cluster_for_mapping %in% keep_clusters]
if (length(keep_cells) < 100) stop("参考核心状态细胞数过少，请检查cluster列。")
obj <- subset(obj, cells = keep_cells)

if ("Assay5" %in% class(obj[["RNA"]]) && exists("JoinLayers", mode = "function")) {
  obj <- tryCatch(
    JoinLayers(obj, assay = "RNA"),
    error = function(e) {
      warning("JoinLayers失败，继续尝试当前RNA层：", e$message)
      obj
    }
  )
}

DefaultAssay(obj) <- "RNA"
logdata <- get_assay_matrix_safe(obj, "data")
if (is.null(logdata) || nrow(logdata) == 0) {
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  logdata <- get_assay_matrix_safe(obj, "data")
}
if (is.null(logdata) || nrow(logdata) == 0) stop("无法取得参考RNA标准化表达矩阵。")

# 统一基因symbol为大写；若大写后重复，则保留平均表达更高的一行。
original_genes <- rownames(logdata)
upper_genes <- toupper(original_genes)
row_means_all <- Matrix::rowMeans(logdata)
ord <- order(row_means_all, decreasing = TRUE)
keep_unique <- !duplicated(upper_genes[ord])
selected_rows <- ord[keep_unique]
logdata <- logdata[selected_rows, , drop = FALSE]
rownames(logdata) <- upper_genes[selected_rows]

state_map <- c(
  C7 = "C7_quiescent",
  C3 = "C3_early_activation",
  C5 = "C5_transition",
  C1 = "C1C2_late_remodeling",
  C2 = "C1C2_late_remodeling",
  C4 = "C4_independent_ECM"
)
state <- unname(state_map[as.character(obj$cluster_for_mapping)])
state_levels <- c(
  "C7_quiescent", "C3_early_activation", "C5_transition",
  "C1C2_late_remodeling", "C4_independent_ECM"
)
state <- factor(state, levels = state_levels)

state_counts <- table(state)
write.csv(
  data.frame(state = names(state_counts), n_reference_cells = as.integer(state_counts)),
  file.path(outdir, "01_reference_state_cell_counts.csv"),
  row.names = FALSE
)

avg_mat <- matrix(NA_real_, nrow = nrow(logdata), ncol = length(state_levels),
                  dimnames = list(rownames(logdata), state_levels))
pct_mat <- avg_mat
for (s in state_levels) {
  idx <- which(state == s)
  avg_mat[, s] <- Matrix::rowMeans(logdata[, idx, drop = FALSE])
  pct_mat[, s] <- Matrix::rowMeans(logdata[, idx, drop = FALSE] > 0)
}

bad_pattern <- grepl("^(MT-|RPL|RPS|MIR|SNOR|HIST|HBA|HBB)", rownames(avg_mat))
classifier_excluded <- rownames(avg_mat) %in% excluded_classifier_genes | bad_pattern

marker_rows <- list()
for (s in state_levels) {
  other_states <- setdiff(state_levels, s)
  specificity_mean <- avg_mat[, s] - rowMeans(avg_mat[, other_states, drop = FALSE])
  specificity_max <- avg_mat[, s] - apply(avg_mat[, other_states, drop = FALSE], 1, max)
  score <- 0.65 * specificity_max + 0.25 * specificity_mean + 0.10 * pct_mat[, s]

  tab <- data.frame(
    gene = rownames(avg_mat),
    state = s,
    avg_logexpr_state = avg_mat[, s],
    pct_state = pct_mat[, s],
    specificity_vs_mean_other = specificity_mean,
    specificity_vs_max_other = specificity_max,
    marker_rank_score = score,
    excluded_from_classifier = classifier_excluded,
    stringsAsFactors = FALSE
  )

  tab <- tab[
    !tab$excluded_from_classifier &
      tab$pct_state >= MIN_PCT &
      tab$avg_logexpr_state >= MIN_AVG_LOGEXPR &
      tab$specificity_vs_max_other > 0,
    , drop = FALSE
  ]
  tab <- tab[order(tab$marker_rank_score, decreasing = TRUE), , drop = FALSE]
  tab <- head(tab, TOP_MARKERS_PER_STATE)

  anchors <- intersect(manual_anchors[[s]], rownames(avg_mat))
  anchors <- setdiff(anchors, tab$gene)
  if (length(anchors) > 0) {
    anchor_tab <- data.frame(
      gene = anchors,
      state = s,
      avg_logexpr_state = avg_mat[anchors, s],
      pct_state = pct_mat[anchors, s],
      specificity_vs_mean_other = avg_mat[anchors, s] - rowMeans(avg_mat[anchors, other_states, drop = FALSE]),
      specificity_vs_max_other = avg_mat[anchors, s] - apply(avg_mat[anchors, other_states, drop = FALSE], 1, max),
      marker_rank_score = NA_real_,
      excluded_from_classifier = FALSE,
      stringsAsFactors = FALSE
    )
    tab <- rbind(tab, anchor_tab)
  }

  tab$marker_source <- ifelse(tab$gene %in% manual_anchors[[s]], "reference_plus_manual_anchor", "reference_specificity")
  marker_rows[[s]] <- tab
}

markers <- do.call(rbind, marker_rows)
markers <- markers[!duplicated(paste(markers$state, markers$gene, sep = "|||")), , drop = FALSE]
write.csv(markers, file.path(outdir, "02_reference_state_classifier_markers.csv"), row.names = FALSE)

classifier_genes <- unique(markers$gene)
centroids <- avg_mat[classifier_genes, state_levels, drop = FALSE]
centroid_df <- data.frame(gene = rownames(centroids), centroids, check.names = FALSE)
write.csv(centroid_df, file.path(outdir, "03_reference_state_centroids.csv"), row.names = FALSE)

availability <- data.frame(
  gene = unique(c(candidate_genes, unlist(manual_anchors))),
  present_in_reference = unique(c(candidate_genes, unlist(manual_anchors))) %in% rownames(avg_mat),
  used_as_classifier_feature = unique(c(candidate_genes, unlist(manual_anchors))) %in% classifier_genes,
  stringsAsFactors = FALSE
)
write.csv(availability, file.path(outdir, "04_reference_target_gene_availability.csv"), row.names = FALSE)

# 参考状态的候选基因平均表达，仅用于后续对照，不进入分类器。
target_present <- intersect(candidate_genes, rownames(avg_mat))
if (length(target_present) > 0) {
  target_ref <- data.frame(gene = target_present, avg_mat[target_present, state_levels, drop = FALSE], check.names = FALSE)
  write.csv(target_ref, file.path(outdir, "05_reference_candidate_state_profile.csv"), row.names = FALSE)
}

summary_lines <- c(
  "Reference state model built from the 53-sample Seurat object.",
  paste0("Reference cluster column: ", cluster_col),
  paste0("Reference cells used: ", ncol(obj)),
  paste0("Classifier genes: ", length(classifier_genes)),
  paste0("States: ", paste(state_levels, collapse = ", ")),
  "C1 and C2 were merged as C1C2_late_remodeling.",
  "C4 was retained as an independent ECM/venous-capillary-like decoy class.",
  "The original nine candidates and SULF1 were excluded from classifier features.",
  "Therefore, later candidate/SULF1 comparisons are not circularly defined by those genes."
)
writeLines(summary_lines, file.path(outdir, "README_reference_model.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_reference_R.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
