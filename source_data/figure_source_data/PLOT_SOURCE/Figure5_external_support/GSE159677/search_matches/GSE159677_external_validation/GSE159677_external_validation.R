options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 30 * 1024^3)
set.seed(20260731)

BASE <- "/root/GSE159677_external_validation"
AGG_DIR <- file.path(BASE, "aggregate_matrix")
OUT_DIR <- file.path(BASE, "results")
REF_EC <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------- package checks ------------------------------------
required <- c("Seurat", "Matrix", "ggplot2", "dplyr", "tidyr", "patchwork", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  cran_missing <- setdiff(missing, c("Seurat", "Matrix"))
  if (length(cran_missing) > 0) {
    try(install.packages(cran_missing, repos = "https://cloud.r-project.org", Ncpus = 4), silent = TRUE)
  }
}
missing2 <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing2) > 0) {
  stop("Missing required R packages: ", paste(missing2, collapse = ", "),
       "\nPlease run this script in the same R 4.3 environment used for the previous Seurat analysis.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
})

# Pure-white scientific figure theme
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
  ggsave(file.path(OUT_DIR, filename), plot = plot, width = width, height = height,
         dpi = dpi, bg = "white", limitsize = FALSE)
}

# Seurat v4/v5 compatibility
join_layers_safe <- function(obj, assay = "RNA") {
  if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = assay), error = function(e) obj)
  } else if (exists("JoinLayers", where = asNamespace("Seurat"), inherits = FALSE)) {
    obj <- tryCatch(Seurat::JoinLayers(obj, assay = assay), error = function(e) obj)
  }
  obj
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

read_10x_any_compression <- function(sample_root, sample_id = "aggregate", prefix_cells = FALSE) {
  matrix_file <- find_one_file(sample_root, c("matrix\\.mtx(\\.gz|\\.bz2|\\.bz)?$"))
  feature_file <- find_one_file(sample_root, c("features\\.tsv(\\.gz|\\.bz2|\\.bz)?$",
                                                "genes\\.tsv(\\.gz|\\.bz2|\\.bz)?$"))
  barcode_file <- find_one_file(sample_root, c("barcodes\\.tsv(\\.gz|\\.bz2|\\.bz)?$"))

  message("Reading ", sample_id)
  message("  matrix:   ", matrix_file)
  message("  features: ", feature_file)
  message("  barcodes: ", barcode_file)

  cm <- open_text_connection(matrix_file)
  on.exit(close(cm), add = TRUE)
  mat <- Matrix::readMM(cm)

  cf <- open_text_connection(feature_file)
  features <- utils::read.delim(cf, header = FALSE, sep = "\t", quote = "", comment.char = "",
                                check.names = FALSE)
  close(cf)

  cb <- open_text_connection(barcode_file)
  barcodes <- utils::read.delim(cb, header = FALSE, sep = "\t", quote = "", comment.char = "",
                                check.names = FALSE)
  close(cb)

  if (nrow(features) != nrow(mat)) {
    stop(sample_id, ": feature rows (", nrow(features), ") != matrix rows (", nrow(mat), ")")
  }
  if (nrow(barcodes) != ncol(mat)) {
    stop(sample_id, ": barcode rows (", nrow(barcodes), ") != matrix columns (", ncol(mat), ")")
  }

  gene_col <- if (ncol(features) >= 2) 2 else 1
  gene_names <- as.character(features[[gene_col]])
  gene_names[is.na(gene_names) | gene_names == ""] <- as.character(features[[1]])[is.na(gene_names) | gene_names == ""]
  gene_names <- make.unique(gene_names)
  raw_barcodes <- as.character(barcodes[[1]])
  cell_names <- if (isTRUE(prefix_cells)) paste(sample_id, raw_barcodes, sep = "_") else raw_barcodes

  mat <- as(mat, "dgCMatrix")
  rownames(mat) <- gene_names
  colnames(mat) <- cell_names
  mat
}

# ------------------------- sample metadata -----------------------------------
meta <- read.csv(file.path(BASE, "sample_metadata.csv"), check.names = FALSE)
meta$condition <- factor(meta$condition, levels = c("PA", "AC"))
write.csv(meta, file.path(OUT_DIR, "00_sample_metadata.csv"), row.names = FALSE)

# ------------------------- read aggregate and split by GEM suffix ------------
agg_mat <- read_10x_any_compression(AGG_DIR, sample_id = "aggregate", prefix_cells = FALSE)

barcode_suffix <- sub("^.*-([0-9]+)$", "\\1", colnames(agg_mat))
if (any(!grepl("^[1-6]$", barcode_suffix))) {
  bad <- unique(barcode_suffix[!grepl("^[1-6]$", barcode_suffix)])
  stop(
    "Unexpected GEM-well suffixes in aggregate barcodes: ",
    paste(head(bad, 20), collapse = ", "),
    "\nExpected suffixes -1 through -6."
  )
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

  if (length(cells) == 0) {
    stop("No cells found for GEM suffix ", suf, " (", sid, ")")
  }

  mat <- agg_mat[, cells, drop = FALSE]
  colnames(mat) <- paste0(sid, "_", colnames(mat))

  obj <- CreateSeuratObject(
    counts = mat,
    project = "GSE159677",
    min.cells = 3,
    min.features = 200
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

# QC violin before downstream integration
qc_meta <- bind_rows(lapply(obj_list, function(x) {
  data.frame(sample_id = x$sample_id, condition = x$condition,
             nFeature_RNA = x$nFeature_RNA, nCount_RNA = x$nCount_RNA,
             percent.mt = x$percent.mt)
}))

p_qc1 <- ggplot(qc_meta, aes(x = sample_id, y = nFeature_RNA, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  labs(x = NULL, y = "Detected genes", title = "GSE159677 QC: genes per cell") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_qc2 <- ggplot(qc_meta, aes(x = sample_id, y = percent.mt, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  labs(x = NULL, y = "Mitochondrial percentage", title = "GSE159677 QC: mitochondrial fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png(p_qc1 / p_qc2, "01_QC_violin_plots.png", 12, 10)

# ------------------------- merge and cluster full dataset --------------------
full <- merge(x = obj_list[[1]], y = obj_list[-1], project = "GSE159677")
full <- join_layers_safe(full, "RNA")
DefaultAssay(full) <- "RNA"
full <- NormalizeData(full, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
full <- FindVariableFeatures(full, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
full <- ScaleData(full, features = VariableFeatures(full), verbose = FALSE)
full <- RunPCA(full, features = VariableFeatures(full), npcs = 30, verbose = FALSE)

reduction_use <- "pca"
if (requireNamespace("harmony", quietly = TRUE)) {
  full <- tryCatch({
    full <- harmony::RunHarmony(full, group.by.vars = "sample_id", reduction = "pca",
                                assay.use = "RNA", verbose = FALSE)
    reduction_use <- "harmony"
    full
  }, error = function(e) {
    message("Harmony failed; using PCA instead: ", conditionMessage(e))
    full
  })
}

full <- RunUMAP(full, reduction = reduction_use, dims = 1:30, seed.use = 20260731, verbose = FALSE)
full <- FindNeighbors(full, reduction = reduction_use, dims = 1:30, verbose = FALSE)
full <- FindClusters(full, resolution = 0.5, random.seed = 20260731, verbose = FALSE)

p_full_condition <- DimPlot(full, reduction = "umap", group.by = "condition", pt.size = 0.15) +
  ggtitle("GSE159677: AC versus PA") + white_theme
p_full_sample <- DimPlot(full, reduction = "umap", group.by = "sample_id", pt.size = 0.12) +
  ggtitle("GSE159677: six paired samples") + white_theme
p_full_cluster <- DimPlot(full, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                          repel = TRUE, pt.size = 0.12) +
  ggtitle("GSE159677: unsupervised clusters") + white_theme
save_png(p_full_condition, "02_full_UMAP_condition.png", 9, 7)
save_png(p_full_sample, "03_full_UMAP_sample.png", 10, 7)
save_png(p_full_cluster, "04_full_UMAP_clusters.png", 10, 8)

# ------------------------- broad cell-lineage scoring ------------------------
lineage_sets <- list(
  EC_core = c("PECAM1", "VWF", "CDH5", "KDR", "ESAM", "EMCN", "ENG", "RAMP2", "PLVAP", "EGFL7", "ERG", "EPAS1"),
  Immune = c("PTPRC", "LST1", "TYROBP", "FCER1G", "LYZ", "CD3D", "CD3E", "MS4A1", "CD79A"),
  SMC_Pericyte = c("ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "CSPG4", "MCAM", "PDGFRB"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PDGFRA", "COL6A1", "C7"),
  Mast = c("TPSAB1", "TPSB2", "KIT", "CPA3"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19")
)
lineage_sets <- lapply(lineage_sets, function(x) intersect(x, rownames(full)))
lineage_sets <- lineage_sets[lengths(lineage_sets) >= 2]

for (nm in names(lineage_sets)) {
  full <- AddModuleScore(full, features = list(lineage_sets[[nm]]), name = paste0(nm, "_score"), seed = 20260731)
  old <- paste0(nm, "_score1")
  full[[nm]] <- full[[old]][, 1]
}

# Cluster-level summaries
cluster_ids <- sort(unique(as.character(full$seurat_clusters)))
cluster_summary <- lapply(cluster_ids, function(cl) {
  cells <- colnames(full)[as.character(full$seurat_clusters) == cl]
  md <- full@meta.data[cells, , drop = FALSE]
  data.frame(
    cluster = cl,
    n_cells = length(cells),
    EC_core = mean(md$EC_core),
    Immune = mean(md$Immune),
    SMC_Pericyte = mean(md$SMC_Pericyte),
    Fibroblast = mean(md$Fibroblast)
  )
}) %>% bind_rows()

norm_data <- get_assay_matrix(full, "RNA", "data")
ec_detection_genes <- intersect(c("PECAM1", "VWF", "CDH5", "KDR", "ESAM", "EMCN", "ENG", "RAMP2", "PLVAP", "EGFL7"), rownames(norm_data))
for (cl in cluster_ids) {
  cells <- colnames(full)[as.character(full$seurat_clusters) == cl]
  det <- if (length(ec_detection_genes) > 0) {
    mean(Matrix::colSums(norm_data[ec_detection_genes, cells, drop = FALSE] > 0) >= 2)
  } else 0
  cluster_summary$frac_cells_with_2plus_EC_markers[cluster_summary$cluster == cl] <- det
}

# Broad, deliberately inclusive EC candidate rule; reference mapping will refine it.
cluster_summary <- cluster_summary %>%
  mutate(
    max_nonEC = pmax(Immune, SMC_Pericyte, Fibroblast),
    EC_candidate = (frac_cells_with_2plus_EC_markers >= 0.10 & EC_core >= max_nonEC - 0.10) |
      (frac_cells_with_2plus_EC_markers >= 0.25)
  )
write.csv(cluster_summary, file.path(OUT_DIR, "02_cluster_lineage_score_summary.csv"), row.names = FALSE)

ec_clusters <- cluster_summary$cluster[cluster_summary$EC_candidate]
if (length(ec_clusters) == 0) stop("No EC candidate clusters were detected. Review 02_cluster_lineage_score_summary.csv")
message("Broad EC candidate clusters: ", paste(ec_clusters, collapse = ", "))

marker_panel_full <- unique(c(
  "PECAM1", "VWF", "CDH5", "KDR", "ESAM", "EMCN", "RAMP2", "PLVAP",
  "PTPRC", "LST1", "TYROBP", "CD3D", "MS4A1",
  "ACTA2", "TAGLN", "MYH11", "RGS5",
  "COL1A1", "DCN", "LUM", "COL3A1"
))
marker_panel_full <- intersect(marker_panel_full, rownames(full))
p_dot_full <- DotPlot(full, features = marker_panel_full, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("Major-lineage markers by cluster") + white_theme
save_png(p_dot_full, "05_full_cluster_marker_dotplot.png", 15, 8)

# ------------------------- EC subsetting and reclustering --------------------
ec_cells <- colnames(full)[as.character(full$seurat_clusters) %in% ec_clusters]
if (length(ec_cells) < 50) stop("Too few EC candidate cells: ", length(ec_cells))
ec <- subset(full, cells = ec_cells)
ec <- join_layers_safe(ec, "RNA")
DefaultAssay(ec) <- "RNA"
ec <- NormalizeData(ec, verbose = FALSE)
ec <- FindVariableFeatures(ec, nfeatures = min(2500, nrow(ec)), verbose = FALSE)
ec <- ScaleData(ec, features = VariableFeatures(ec), verbose = FALSE)
ec <- RunPCA(ec, npcs = 30, verbose = FALSE)

red_ec <- "pca"
if (requireNamespace("harmony", quietly = TRUE)) {
  ec <- tryCatch({
    ec <- harmony::RunHarmony(ec, group.by.vars = "sample_id", reduction = "pca",
                              assay.use = "RNA", verbose = FALSE)
    red_ec <- "harmony"
    ec
  }, error = function(e) {
    message("EC Harmony failed; using PCA: ", conditionMessage(e))
    ec
  })
}
ec <- RunUMAP(ec, reduction = red_ec, dims = 1:min(30, ncol(Embeddings(ec, red_ec))),
              seed.use = 20260731, verbose = FALSE)
ec <- FindNeighbors(ec, reduction = red_ec, dims = 1:min(30, ncol(Embeddings(ec, red_ec))), verbose = FALSE)
ec <- FindClusters(ec, resolution = 0.4, random.seed = 20260731, verbose = FALSE)

p_ec_condition <- DimPlot(ec, group.by = "condition", pt.size = 0.35) +
  ggtitle("External EC candidates: AC versus PA") + white_theme
p_ec_sample <- DimPlot(ec, group.by = "sample_id", pt.size = 0.30) +
  ggtitle("External EC candidates by sample") + white_theme
p_ec_cluster <- DimPlot(ec, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.30) +
  ggtitle("External EC reclustering") + white_theme
save_png(p_ec_condition, "06_EC_UMAP_condition.png", 9, 7)
save_png(p_ec_sample, "07_EC_UMAP_sample.png", 10, 7)
save_png(p_ec_cluster, "08_EC_UMAP_clusters.png", 9, 7)

# EC counts by sample
cell_counts <- as.data.frame(table(ec$sample_id, ec$condition), stringsAsFactors = FALSE)
colnames(cell_counts) <- c("sample_id", "condition", "n_EC_candidates")
cell_counts <- cell_counts[cell_counts$n_EC_candidates > 0, ]
write.csv(cell_counts, file.path(OUT_DIR, "03_EC_candidate_cell_counts.csv"), row.names = FALSE)
p_ec_counts <- ggplot(cell_counts, aes(x = sample_id, y = n_EC_candidates, fill = condition)) +
  geom_col(width = 0.75) +
  labs(x = NULL, y = "EC candidate cells", title = "EC candidate cell counts by sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png(p_ec_counts, "09_EC_cell_counts_by_sample.png", 9, 6)

# ------------------------- discovery reference preparation ------------------
# Fixed biological anchors; reference-derived markers will be added when possible.
anchor_sets <- list(
  C7_root_quiescent = c("CXCL12", "EFNB2", "GJA4", "GJA5", "SOX17", "HEY1", "KLF2", "KLF4", "ERG"),
  C3_shared_activation_hub = c("ANGPT2", "CCL14", "PGF", "ESM1", "KDR", "EDNRB", "PLVAP"),
  C1_C5_inflammatory_adhesion = c("VCAM1", "SELE", "ICAM1", "CXCL2", "CXCL3", "NFKBIA", "IFITM3", "IFIT1", "ISG15", "STAT1"),
  C2_vascular_remodeling = c("ITLN1", "SULF1", "MGP", "CLU", "ELN", "COL8A1", "COL8A2", "EFEMP1", "MFAP4"),
  C4_ECM_EndMT_like = c("POSTN", "IGFBP5", "SPARC", "COL6A2", "COL1A1", "COL1A2", "COL3A1", "TIMP1", "MAFB", "NR2F2")
)

if (!file.exists(REF_EC)) {
  ref_hits <- tryCatch(
    list.files("/root", pattern = "slingshot.*scored\\.rds$|arterial_dysfunction.*\\.rds$",
               recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
    error = function(e) character(0)
  )
  if (length(ref_hits) > 0) {
    REF_EC <- ref_hits[1]
    message("Fixed reference path was absent; automatically found: ", REF_EC)
  }
}
ref_available <- file.exists(REF_EC)
ref_markers <- data.frame()
state_transfer_success <- FALSE

normalize_cluster_label <- function(x) {
  x <- as.character(x)
  x <- gsub("^cluster[_ ]*", "", x, ignore.case = TRUE)
  x <- gsub("^C", "", x, ignore.case = TRUE)
  x
}

find_reference_cluster_labels <- function(ref) {
  wanted <- c("1", "2", "3", "4", "5", "7")
  idv <- normalize_cluster_label(Idents(ref))
  if (all(wanted %in% unique(idv))) return(idv)
  for (cn in colnames(ref@meta.data)) {
    vv <- normalize_cluster_label(ref@meta.data[[cn]])
    if (all(wanted %in% unique(vv))) {
      message("Using reference cluster metadata column: ", cn)
      return(vv)
    }
  }
  stop("Could not identify reference EC cluster labels containing 1,2,3,4,5,7")
}

if (ref_available) {
  message("Loading discovery EC reference: ", REF_EC)
  ref <- readRDS(REF_EC)
  if (!"RNA" %in% Assays(ref)) stop("Reference object has no RNA assay")
  ref <- join_layers_safe(ref, "RNA")
  DefaultAssay(ref) <- "RNA"
  ref$reference_cluster_num <- find_reference_cluster_labels(ref)
  ref <- subset(ref, subset = reference_cluster_num %in% c("1", "2", "3", "4", "5", "7"))
  ref$reference_state <- paste0("C", ref$reference_cluster_num)
  Idents(ref) <- ref$reference_state

  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, nfeatures = min(3000, nrow(ref)), verbose = FALSE)
  ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
  ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 30, verbose = FALSE)

  # Reference-specific markers for fixed state programs
  ref_markers <- tryCatch({
    FindAllMarkers(ref, assay = "RNA", only.pos = TRUE, min.pct = 0.15,
                   logfc.threshold = 0.20, return.thresh = 0.05, verbose = FALSE)
  }, error = function(e) {
    message("Reference marker discovery failed: ", conditionMessage(e))
    data.frame()
  })
  if (nrow(ref_markers) > 0) {
    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(ref_markers))[1]
    if (!is.na(fc_col)) {
      ref_markers <- ref_markers %>%
        filter(!grepl("^(MT-|RPL|RPS|HBA|HBB)", gene)) %>%
        arrange(cluster, desc(.data[[fc_col]]))
    }
    write.csv(ref_markers, file.path(OUT_DIR, "04_reference_state_markers_all.csv"), row.names = FALSE)
  }

  # State transfer to external EC candidates
  common_features <- intersect(VariableFeatures(ref), rownames(ec))
  if (length(common_features) >= 500) {
    anchors <- tryCatch(
      FindTransferAnchors(
        reference = ref,
        query = ec,
        normalization.method = "LogNormalize",
        reference.reduction = "pca",
        features = common_features,
        dims = 1:30,
        verbose = FALSE
      ),
      error = function(e) {
        message("FindTransferAnchors failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(anchors)) {
      pred <- tryCatch(
        TransferData(anchorset = anchors, refdata = ref$reference_state,
                     dims = 1:30, verbose = FALSE),
        error = function(e) {
          message("TransferData failed: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(pred)) {
        ec <- AddMetaData(ec, pred)
        ec$mapped_state <- ifelse(ec$prediction.score.max >= 0.40,
                                  as.character(ec$predicted.id), "Uncertain")
        state_transfer_success <- TRUE
      }
    }
  }
}

# Add top reference markers to each fixed program, while retaining biological anchors.
if (nrow(ref_markers) > 0) {
  top_ref <- ref_markers %>%
    group_by(cluster) %>%
    slice_head(n = 20) %>%
    ungroup()
  state_to_program <- c(
    C7 = "C7_root_quiescent",
    C3 = "C3_shared_activation_hub",
    C1 = "C1_C5_inflammatory_adhesion",
    C5 = "C1_C5_inflammatory_adhesion",
    C2 = "C2_vascular_remodeling",
    C4 = "C4_ECM_EndMT_like"
  )
  for (st in names(state_to_program)) {
    genes_st <- unique(top_ref$gene[as.character(top_ref$cluster) == st])
    pg <- state_to_program[[st]]
    anchor_sets[[pg]] <- unique(c(anchor_sets[[pg]], genes_st))
  }
}

# Keep available genes and save final signatures
anchor_sets <- lapply(anchor_sets, function(x) unique(intersect(x, rownames(ec))))
sig_table <- bind_rows(lapply(names(anchor_sets), function(nm) {
  data.frame(program = nm, gene = anchor_sets[[nm]])
}))
write.csv(sig_table, file.path(OUT_DIR, "05_fixed_validation_gene_signatures.csv"), row.names = FALSE)

bad_sig <- names(anchor_sets)[lengths(anchor_sets) < 3]
if (length(bad_sig) > 0) warning("Signatures with <3 available genes: ", paste(bad_sig, collapse = ", "))

# Module scores
for (nm in names(anchor_sets)) {
  if (length(anchor_sets[[nm]]) >= 2) {
    ec <- AddModuleScore(ec, features = list(anchor_sets[[nm]]), name = paste0(nm, "_score"), seed = 20260731)
    ec[[nm]] <- ec[[paste0(nm, "_score1")]][, 1]
  } else {
    ec[[nm]] <- NA_real_
  }
}

if (!state_transfer_success) {
  # Fallback: highest fixed program, explicitly marked as program-based rather than label-transfer-based.
  score_mat <- ec@meta.data[, names(anchor_sets), drop = FALSE]
  best_program <- apply(score_mat, 1, function(z) {
    if (all(is.na(z))) return("Uncertain")
    names(z)[which.max(z)]
  })
  fallback_map <- c(
    C7_root_quiescent = "C7",
    C3_shared_activation_hub = "C3",
    C1_C5_inflammatory_adhesion = "C1_C5",
    C2_vascular_remodeling = "C2",
    C4_ECM_EndMT_like = "C4"
  )
  ec$mapped_state <- unname(fallback_map[best_program])
  ec$mapped_state[is.na(ec$mapped_state)] <- "Uncertain"
  ec$prediction.score.max <- NA_real_
}

# Save state mapping diagnostics
map_diag <- ec@meta.data %>%
  count(sample_id, patient_id, condition, mapped_state, name = "n_cells") %>%
  group_by(sample_id) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(map_diag, file.path(OUT_DIR, "06_mapped_state_counts_and_proportions.csv"), row.names = FALSE)

p_ec_state <- DimPlot(ec, group.by = "mapped_state", pt.size = 0.35, label = FALSE) +
  ggtitle("External EC: mapped discovery states") + white_theme
if (state_transfer_success) {
  p_ec_conf <- FeaturePlot(ec, features = "prediction.score.max", pt.size = 0.35) +
    ggtitle("Reference-mapping confidence") + white_theme
} else {
  p_ec_conf <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Reference label transfer was unavailable;\nstate assignment used fixed program scores.", size = 5) +
    xlim(-1, 1) + ylim(-1, 1) + labs(title = "Reference-mapping confidence") +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
}
save_png(p_ec_state, "10_EC_UMAP_mapped_state.png", 10, 7)
save_png(p_ec_conf, "11_EC_UMAP_mapping_confidence.png", 9, 7)

# Marker dotplot
state_marker_panel <- unique(c(
  "PECAM1", "VWF", "CDH5", "KDR", "RAMP2",
  "CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1",
  "ANGPT2", "CCL14", "PGF",
  "VCAM1", "SELE", "ICAM1", "ISG15", "STAT1",
  "ITLN1", "SULF1", "MGP", "CLU", "ELN", "COL8A1",
  "POSTN", "IGFBP5", "SPARC", "COL6A2", "MAFB", "NR2F2"
))
state_marker_panel <- intersect(state_marker_panel, rownames(ec))
Idents(ec) <- ec$mapped_state
p_state_dot <- DotPlot(ec, features = state_marker_panel, group.by = "mapped_state") +
  RotatedAxis() + ggtitle("Fixed discovery markers in GSE159677 EC") + white_theme
save_png(p_state_dot, "12_EC_state_marker_dotplot.png", 17, 8)

# ------------------------- sample-level program validation ------------------
program_names <- names(anchor_sets)
sample_program <- ec@meta.data %>%
  group_by(sample_id, patient_id, condition) %>%
  summarise(
    n_EC = n(),
    across(all_of(program_names), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
write.csv(sample_program, file.path(OUT_DIR, "07_sample_level_program_scores.csv"), row.names = FALSE)

program_long <- sample_program %>%
  pivot_longer(cols = all_of(program_names), names_to = "program", values_to = "mean_score")

# Paired plot: points from same patient connected
p_program_pair <- ggplot(program_long,
                         aes(x = condition, y = mean_score, group = patient_id, colour = patient_id)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ program, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Mean EC program score",
       title = "Paired external validation: PA versus AC") +
  white_theme
save_png(p_program_pair, "13_paired_program_scores_PA_vs_AC.png", 12, 11)

# Heatmap of sample-level scores
heat_mat <- as.data.frame(sample_program)
rownames(heat_mat) <- heat_mat$sample_id
heat_mat <- as.matrix(heat_mat[, program_names, drop = FALSE])
heat_z <- t(scale(t(heat_mat)))
heat_z[!is.finite(heat_z)] <- 0
png(file.path(OUT_DIR, "14_sample_program_score_heatmap.png"), width = 3000, height = 2200, res = 300, bg = "white")
pheatmap::pheatmap(heat_z, cluster_rows = TRUE, cluster_cols = FALSE,
                   main = "GSE159677 EC program scores (row z-score)", border_color = NA)
dev.off()

# Paired statistics helper
paired_stats <- function(df, value_col, feature_col, expected_map = NULL) {
  features <- unique(df[[feature_col]])
  out <- lapply(features, function(feat) {
    dd <- df[df[[feature_col]] == feat, c("patient_id", "condition", value_col), drop = FALSE]
    wide <- tidyr::pivot_wider(dd, names_from = condition, values_from = all_of(value_col))
    wide <- wide[complete.cases(wide[, c("PA", "AC")]), ]
    diff <- wide$AC - wide$PA
    p <- if (length(diff) >= 2) {
      tryCatch(wilcox.test(wide$AC, wide$PA, paired = TRUE, exact = TRUE)$p.value,
               error = function(e) NA_real_)
    } else NA_real_
    dz <- if (length(diff) >= 2 && sd(diff) > 0) mean(diff) / sd(diff) else NA_real_
    expected <- if (!is.null(expected_map) && feat %in% names(expected_map)) expected_map[[feat]] else NA_character_
    n_expected <- if (is.na(expected)) NA_integer_ else if (expected == "AC>PA") sum(diff > 0) else sum(diff < 0)
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
      P1_AC_minus_PA = diff[match("P1", wide$patient_id)],
      P2_AC_minus_PA = diff[match("P2", wide$patient_id)],
      P3_AC_minus_PA = diff[match("P3", wide$patient_id)]
    )
  })
  bind_rows(out)
}

expected_program <- c(
  C7_root_quiescent = "PA>AC",
  C3_shared_activation_hub = "AC>PA",
  C1_C5_inflammatory_adhesion = "AC>PA",
  C2_vascular_remodeling = "AC>PA",
  C4_ECM_EndMT_like = "AC>PA"
)
program_stats <- paired_stats(program_long, "mean_score", "program", expected_program)
write.csv(program_stats, file.path(OUT_DIR, "08_paired_program_statistics.csv"), row.names = FALSE)

# ------------------------- mapped-state proportions -------------------------
state_prop <- map_diag
p_state_pair <- ggplot(state_prop,
                       aes(x = condition, y = proportion, group = patient_id, colour = patient_id)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ mapped_state, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Proportion among EC candidates",
       title = "Paired mapped-state proportions") +
  white_theme
save_png(p_state_pair, "15_paired_mapped_state_proportions.png", 13, 9)

state_expected <- c(C7 = "PA>AC", C3 = "AC>PA", C1 = "AC>PA", C5 = "AC>PA",
                    C1_C5 = "AC>PA", C2 = "AC>PA", C4 = "AC>PA", Uncertain = NA_character_)
state_stats <- paired_stats(state_prop, "proportion", "mapped_state", state_expected)
write.csv(state_stats, file.path(OUT_DIR, "09_paired_state_proportion_statistics.csv"), row.names = FALSE)

# ------------------------- key-gene sample summaries ------------------------
key_genes <- intersect(c(
  "CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1",
  "ANGPT2", "CCL14", "PGF",
  "VCAM1", "SELE", "ICAM1",
  "ITLN1", "SULF1", "MGP", "CLU", "ELN", "COL8A1",
  "POSTN", "IGFBP5", "SPARC", "COL6A2", "MAFB", "NR2F2"
), rownames(ec))

edata <- get_assay_matrix(ec, "RNA", "data")
key_gene_rows <- lapply(unique(ec$sample_id), function(sid) {
  cells <- colnames(ec)[ec$sample_id == sid]
  md1 <- ec@meta.data[cells[1], c("sample_id", "patient_id", "condition"), drop = FALSE]
  vals <- Matrix::rowMeans(edata[key_genes, cells, drop = FALSE])
  data.frame(sample_id = md1$sample_id, patient_id = md1$patient_id,
             condition = md1$condition, gene = names(vals), mean_log_normalized_expression = as.numeric(vals))
})
key_gene_sample <- bind_rows(key_gene_rows)
write.csv(key_gene_sample, file.path(OUT_DIR, "10_key_gene_sample_mean_expression.csv"), row.names = FALSE)

p_key_gene <- ggplot(key_gene_sample,
                     aes(x = condition, y = mean_log_normalized_expression,
                         group = patient_id, colour = patient_id)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = "Mean log-normalized expression",
       title = "Key genes: paired PA versus AC") +
  white_theme
save_png(p_key_gene, "16_key_gene_paired_expression.png", 16, 13)

key_expected <- setNames(rep("AC>PA", length(key_genes)), key_genes)
key_expected[intersect(names(key_expected), c("CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1"))] <- "PA>AC"
key_stats <- paired_stats(key_gene_sample, "mean_log_normalized_expression", "gene", key_expected)
write.csv(key_stats, file.path(OUT_DIR, "11_key_gene_paired_statistics.csv"), row.names = FALSE)

# ------------------------- pseudobulk counts and edgeR -----------------------
counts <- get_assay_matrix(ec, "RNA", "counts")
samples <- unique(ec$sample_id)
pb_counts <- sapply(samples, function(sid) {
  Matrix::rowSums(counts[, ec$sample_id == sid, drop = FALSE])
})
if (!is.matrix(pb_counts)) pb_counts <- as.matrix(pb_counts)
colnames(pb_counts) <- samples
write.csv(data.frame(gene = rownames(pb_counts), pb_counts, check.names = FALSE),
          file.path(OUT_DIR, "12_EC_pseudobulk_raw_counts.csv"), row.names = FALSE)

# edgeR is optional but scientifically preferable for paired pseudobulk effect estimates.
if (!requireNamespace("edgeR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    try(install.packages("BiocManager", repos = "https://cloud.r-project.org"), silent = TRUE)
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    try(BiocManager::install("edgeR", ask = FALSE, update = FALSE), silent = TRUE)
  }
}

if (requireNamespace("edgeR", quietly = TRUE)) {
  edgeR_result <- tryCatch({
    suppressPackageStartupMessages(library(edgeR))
    pb_meta <- meta[match(colnames(pb_counts), meta$sample_id), ]
    pb_meta$condition <- factor(pb_meta$condition, levels = c("PA", "AC"))
    pb_meta$patient_id <- factor(pb_meta$patient_id)
    design <- model.matrix(~ patient_id + condition, data = pb_meta)
    if (qr(design)$rank < ncol(design)) stop("Paired design matrix is not full rank")
    y <- edgeR::DGEList(counts = pb_counts)
    keep_gene <- edgeR::filterByExpr(y, design = design, min.count = 5)
    y <- y[keep_gene, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    coef_name <- grep("conditionAC", colnames(design), value = TRUE)
    if (length(coef_name) != 1) stop("Could not identify conditionAC coefficient")
    qlf <- edgeR::glmQLFTest(fit, coef = which(colnames(design) == coef_name))
    tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
    tab$gene <- rownames(tab)
    tab <- tab[, c("gene", setdiff(colnames(tab), "gene"))]
    write.csv(tab, file.path(OUT_DIR, "13_EC_pseudobulk_edgeR_paired_AC_vs_PA_all_genes.csv"), row.names = FALSE)
    write.csv(tab[tab$gene %in% key_genes, ],
              file.path(OUT_DIR, "14_EC_pseudobulk_edgeR_key_genes.csv"), row.names = FALSE)
    TRUE
  }, error = function(e) {
    writeLines(paste("edgeR paired analysis failed:", conditionMessage(e)),
               file.path(OUT_DIR, "13_edgeR_NOT_RUN.txt"))
    FALSE
  })
} else {
  writeLines("edgeR was not available; raw pseudobulk counts were generated, but formal edgeR testing was skipped.",
             file.path(OUT_DIR, "13_edgeR_NOT_RUN.txt"))
}

# ------------------------- validation summary table -------------------------
summary_program <- program_stats %>%
  transmute(
    validation_item = feature,
    expected_direction,
    mean_AC_minus_PA,
    n_pairs_in_expected_direction,
    fraction_in_expected_direction,
    evidence_grade = case_when(
      n_pairs_in_expected_direction == 3 ~ "directionally_consistent_3_of_3",
      n_pairs_in_expected_direction == 2 ~ "partially_consistent_2_of_3",
      n_pairs_in_expected_direction <= 1 ~ "not_consistent",
      TRUE ~ "not_assessed"
    )
  )
write.csv(summary_program, file.path(OUT_DIR, "15_primary_validation_summary.csv"), row.names = FALSE)

# Save objects and session information
saveRDS(full, file.path(OUT_DIR, "GSE159677_full_QC_clustered.rds"), compress = FALSE)
saveRDS(ec, file.path(OUT_DIR, "GSE159677_EC_external_validation.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# Human-readable run summary
summary_text <- c(
  "GSE159677 independent external validation completed.",
  paste0("Full cells after QC: ", ncol(full)),
  paste0("Broad EC candidate cells: ", ncol(ec)),
  paste0("EC candidate clusters from full atlas: ", paste(ec_clusters, collapse = ", ")),
  paste0("Discovery reference available: ", ref_available),
  paste0("State label transfer completed: ", state_transfer_success),
  "Primary interpretation should emphasize paired direction consistency and effect size.",
  "With only three patient pairs, Wilcoxon P values have very limited resolution and should not be treated as the sole reliability criterion.",
  "PA is proximal adjacent diseased artery, not a truly healthy control."
)
writeLines(summary_text, file.path(OUT_DIR, "README_results.txt"))
cat(paste(summary_text, collapse = "\n"), "\n")
