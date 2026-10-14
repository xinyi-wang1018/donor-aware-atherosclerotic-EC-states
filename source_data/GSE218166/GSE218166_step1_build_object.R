suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(harmony)
})

set.seed(20260812)

# ============================================================
# 0. Paths
# ============================================================
BASE <- "/root/GSE218166"
OUT  <- file.path(BASE, "step1_build_and_broad_cluster")
CACHE <- file.path(OUT, "cache")
PLOT <- file.path(OUT, "plots")
TAB  <- file.path(OUT, "tables")
LOG  <- file.path(OUT, "logs")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 200 * 1024^3)

nthread <- min(32, parallel::detectCores())
data.table::setDTthreads(nthread)

cat("============================================\n")
cat("GSE218166 STEP1 START\n")
cat("Threads:", nthread, "\n")
cat("============================================\n\n")


# ============================================================
# 1. Helper: read count matrix
#
# GEO matrix format:
# header = cell1 cell2 cell3 ...
# row    = gene count count count ...
#
# Header lacks the first "gene" field.
# We prepend "gene<TAB>" during streaming.
# ============================================================
read_count_matrix <- function(file, prefix, cache_file) {

  if (file.exists(cache_file)) {
    cat("\nLoading cached sparse matrix:\n", cache_file, "\n")
    return(readRDS(cache_file))
  }

  cat("\n============================================\n")
  cat("Reading:", basename(file), "\n")
  cat("============================================\n")

  cmd <- sprintf(
    "bash -c \"printf 'gene\\t'; zcat '%s'\"",
    file
  )

  dt <- fread(
    cmd = cmd,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    showProgress = TRUE,
    nThread = nthread
  )

  cat("Dense table dimensions:", nrow(dt), "x", ncol(dt), "\n")

  gene <- as.character(dt[[1]])

  if (anyNA(gene) || any(gene == "")) {
    stop("Missing gene names detected in ", file)
  }

  ndups <- sum(duplicated(gene))
  cat("Duplicated gene names:", ndups, "\n")

  # Keep first occurrence unchanged; suffix duplicates if present.
  # This avoids Seurat rowname failure while retaining all rows.
  gene_unique <- make.unique(gene)

  cell_raw <- names(dt)[-1]

  cat("Number of cells:", length(cell_raw), "\n")

  # Dense -> matrix
  dense_mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(dense_mat) <- "numeric"

  rownames(dense_mat) <- gene_unique
  colnames(dense_mat) <- paste0(prefix, "__", cell_raw)

  rm(dt)
  gc()

  cat("Converting to sparse matrix...\n")

  sp <- Matrix(
    dense_mat,
    sparse = TRUE
  )

  rm(dense_mat)
  gc()

  sp <- as(sp, "dgCMatrix")

  cat("Sparse dimensions:", nrow(sp), "x", ncol(sp), "\n")
  cat("Non-zero entries:", length(sp@x), "\n")

  saveRDS(
    sp,
    cache_file,
    compress = FALSE
  )

  cat("Saved cache:", cache_file, "\n")

  return(sp)
}


# ============================================================
# 2. Helper: read metadata
# ============================================================
read_metadata <- function(file, prefix) {

  cat("\nReading metadata:", basename(file), "\n")

  md <- fread(
    cmd = paste("zcat", shQuote(file)),
    sep = ",",
    header = TRUE,
    check.names = FALSE
  )

  # First unnamed column = cell ID
  setnames(md, 1, "cell_id_raw")

  # Preserve published QC values using new names
  if ("nCount_RNA" %in% names(md)) {
    setnames(md, "nCount_RNA", "nCount_published")
  }

  if ("nFeature_RNA" %in% names(md)) {
    setnames(md, "nFeature_RNA", "nFeature_published")
  }

  md[, dataset_part := prefix]
  md[, cell_id := paste0(prefix, "__", cell_id_raw)]

  return(md)
}


# ============================================================
# 3. Input files
# ============================================================
FL_COUNT <- file.path(
  BASE,
  "GSM6735811_GeneMatrix_counts.csv.gz"
)

FL_META <- file.path(
  BASE,
  "GSM6735811_CellInfo.csv.gz"
)

DE3_COUNT <- file.path(
  BASE,
  "GSM6753958_sepDE3_GeneMatrix_counts.csv.gz"
)

DE3_META <- file.path(
  BASE,
  "GSM6753958_sepDE3_CellInfo.csv.gz"
)


# ============================================================
# 4. Read / cache sparse matrices
# ============================================================
fl <- read_count_matrix(
  FL_COUNT,
  prefix = "FL",
  cache_file = file.path(CACHE, "FL_counts_sparse.rds")
)

de3 <- read_count_matrix(
  DE3_COUNT,
  prefix = "3DE",
  cache_file = file.path(CACHE, "3DE_counts_sparse.rds")
)

cat("\nChecking gene order...\n")

if (!identical(rownames(fl), rownames(de3))) {
  stop(
    "FL and 3DE gene rows are not identical. ",
    "Do NOT continue automatically."
  )
}

cat("Gene rows identical: YES\n")


# ============================================================
# 5. Combine matrices
# ============================================================
cat("\nCombining FL + 3DE matrices...\n")

counts <- cbind(fl, de3)

rm(fl, de3)
gc()

# Remove genes with absolutely no counts
keep_gene <- Matrix::rowSums(counts) > 0
counts <- counts[keep_gene, , drop = FALSE]

cat(
  "Combined matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n"
)


# ============================================================
# 6. Metadata
# ============================================================
md_fl  <- read_metadata(FL_META, "FL")
md_de3 <- read_metadata(DE3_META, "3DE")

md <- rbindlist(
  list(md_fl, md_de3),
  fill = TRUE,
  use.names = TRUE
)

rm(md_fl, md_de3)
gc()

cat("Metadata rows:", nrow(md), "\n")
cat("Matrix cells:", ncol(counts), "\n")

if (anyDuplicated(md$cell_id)) {
  stop("Duplicated metadata cell IDs detected.")
}

missing_meta <- setdiff(colnames(counts), md$cell_id)
missing_mat  <- setdiff(md$cell_id, colnames(counts))

cat("Cells missing metadata:", length(missing_meta), "\n")
cat("Metadata missing matrix:", length(missing_mat), "\n")

if (
  length(missing_meta) > 0 ||
  length(missing_mat) > 0
) {
  stop("Counts / metadata cell IDs do not perfectly match.")
}

# Match metadata exactly to matrix column order
md <- md[match(colnames(counts), cell_id)]

stopifnot(
  identical(md$cell_id, colnames(counts))
)

# Standardized variables
md[, disease_stage := factor(
  Level_nature2016,
  levels = c("Normal", "Early", "Advanced"),
  ordered = TRUE
)]

md[, segment_id := paste(
  Individual,
  position,
  sep = "_"
)]

md[, donor_id := Individual]
md[, coronary_position := position]
md[, lesion_process := Process_nature2016]

md_df <- as.data.frame(md)
rownames(md_df) <- md_df$cell_id


# ============================================================
# 7. Create Seurat object
#
# Important:
# These are already published processed cells.
# Step1 does NOT apply an aggressive second QC filter.
# ============================================================
cat("\nCreating Seurat object...\n")

obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE218166",
  meta.data = md_df,
  min.cells = 0,
  min.features = 0
)

rm(counts, md, md_df)
gc()

obj[["percent.mt"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

cat(
  "Seurat object:",
  nrow(obj),
  "genes x",
  ncol(obj),
  "cells\n"
)


# ============================================================
# 8. Save metadata / composition checks
# ============================================================
meta <- obj@meta.data

# Cell numbers by stage / technology
tab_stage_tech <- as.data.table(
  table(
    tech = meta$tech,
    stage = meta$disease_stage
  )
)

fwrite(
  tab_stage_tech,
  file.path(TAB, "01_cells_by_tech_stage.csv")
)

# Unique segment information
segment_info <- unique(
  meta[
    ,
    c(
      "segment_id",
      "donor_id",
      "coronary_position",
      "lesion_process",
      "disease_stage",
      "tech"
    )
  ]
)

segment_info <- segment_info[
  order(
    segment_info$tech,
    segment_info$donor_id,
    segment_info$coronary_position
  ),
]

write.csv(
  segment_info,
  file.path(TAB, "02_segment_information.csv"),
  row.names = FALSE
)

# QC summary
qc_dt <- as.data.table(meta)

qc_summary <- qc_dt[
  ,
  .(
    n_cells = .N,

    nCount_median =
      as.numeric(median(nCount_RNA)),

    nFeature_median =
      as.numeric(median(nFeature_RNA)),

    percent_mt_median =
      median(percent.mt),

    percent_mt_Q95 =
      as.numeric(quantile(percent.mt, 0.95, names = FALSE)),

    nFeature_Q05 =
      as.numeric(quantile(nFeature_RNA, 0.05, names = FALSE)),

    nFeature_Q95 =
      as.numeric(quantile(nFeature_RNA, 0.95, names = FALSE))
  ),
  by = .(
    tech,
    disease_stage
  )
]

fwrite(
  qc_summary,
  file.path(TAB, "03_QC_summary_by_stage_tech.csv")
)


# ============================================================
# 9. QC plots
# ============================================================
p_qc1 <- VlnPlot(
  obj,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  ),
  group.by = "tech",
  pt.size = 0,
  ncol = 3
) &
  theme_bw()

ggsave(
  file.path(PLOT, "01_QC_by_tech.png"),
  p_qc1,
  width = 15,
  height = 5,
  dpi = 300,
  bg = "white"
)

p_qc2 <- VlnPlot(
  obj,
  features = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  ),
  group.by = "disease_stage",
  pt.size = 0,
  ncol = 3
) &
  theme_bw()

ggsave(
  file.path(PLOT, "02_QC_by_stage.png"),
  p_qc2,
  width = 15,
  height = 5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 10. Standard RNA workflow
# ============================================================
cat("\nNormalizeData...\n")

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("FindVariableFeatures...\n")

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

cat("ScaleData...\n")

obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  verbose = FALSE
)

cat("RunPCA...\n")

obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 50,
  verbose = FALSE
)


# ============================================================
# 11. Raw PCA UMAP
# Used to inspect technology effects BEFORE correction
# ============================================================
cat("Run raw PCA UMAP...\n")

obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.pca",
  reduction.key = "UMAPPCA_",
  seed.use = 20260812,
  verbose = FALSE
)

p_raw_tech <- DimPlot(
  obj,
  reduction = "umap.pca",
  group.by = "tech"
) +
  theme_bw()

ggsave(
  file.path(PLOT, "03_UMAP_rawPCA_by_tech.png"),
  p_raw_tech,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

p_raw_stage <- DimPlot(
  obj,
  reduction = "umap.pca",
  group.by = "disease_stage"
) +
  theme_bw()

ggsave(
  file.path(PLOT, "04_UMAP_rawPCA_by_stage.png"),
  p_raw_stage,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 12. Harmony
#
# ONLY for:
#   broad clustering / visualization / cell-type identification
#
# NOT for:
#   stage DEG
#   pseudobulk
#   Normal-Early-Advanced biological inference
# ============================================================
cat("\nRunning Harmony on tech...\n")

obj <- harmony::RunHarmony(
  obj,
  group.by.vars = "tech",
  reduction = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  verbose = TRUE
)

cat("Harmony UMAP...\n")

obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap.harmony",
  reduction.key = "UMAPH_",
  seed.use = 20260812,
  verbose = FALSE
)

obj <- FindNeighbors(
  obj,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  resolution = 0.6,
  random.seed = 20260812,
  verbose = FALSE
)


# ============================================================
# 13. UMAP plots
# ============================================================
plots <- list(

  cluster = DimPlot(
    obj,
    reduction = "umap.harmony",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) +
    theme_bw(),

  tech = DimPlot(
    obj,
    reduction = "umap.harmony",
    group.by = "tech"
  ) +
    theme_bw(),

  stage = DimPlot(
    obj,
    reduction = "umap.harmony",
    group.by = "disease_stage"
  ) +
    theme_bw(),

  donor = DimPlot(
    obj,
    reduction = "umap.harmony",
    group.by = "donor_id"
  ) +
    theme_bw(),

  position = DimPlot(
    obj,
    reduction = "umap.harmony",
    group.by = "coronary_position"
  ) +
    theme_bw()
)

ggsave(
  file.path(PLOT, "05_UMAP_Harmony_clusters.png"),
  plots$cluster,
  width = 9,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(PLOT, "06_UMAP_Harmony_tech.png"),
  plots$tech,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(PLOT, "07_UMAP_Harmony_stage.png"),
  plots$stage,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(PLOT, "08_UMAP_Harmony_donor.png"),
  plots$donor,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(PLOT, "09_UMAP_Harmony_coronary_position.png"),
  plots$position,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 14. Cluster composition tables
# ============================================================
cluster_count <- as.data.table(
  table(
    cluster = obj$seurat_clusters
  )
)

fwrite(
  cluster_count,
  file.path(TAB, "04_cluster_cell_counts.csv")
)

cluster_tech <- as.data.table(
  table(
    cluster = obj$seurat_clusters,
    tech = obj$tech
  )
)

fwrite(
  cluster_tech,
  file.path(TAB, "05_cluster_by_tech.csv")
)

cluster_stage <- as.data.table(
  table(
    cluster = obj$seurat_clusters,
    stage = obj$disease_stage
  )
)

fwrite(
  cluster_stage,
  file.path(TAB, "06_cluster_by_stage.csv")
)

cluster_position <- as.data.table(
  table(
    cluster = obj$seurat_clusters,
    position = obj$coronary_position
  )
)

fwrite(
  cluster_position,
  file.path(TAB, "07_cluster_by_coronary_position.csv")
)


# ============================================================
# 15. Broad marker panel
#
# IMPORTANT:
# We do NOT automatically delete fibroblast-like genes from EC.
# Diseased EC may acquire remodeling / mesenchymal programs.
# ============================================================
marker_sets <- list(

  Endothelial = c(
    "PECAM1",
    "VWF",
    "KDR",
    "CLDN5",
    "ESAM",
    "ENG",
    "RAMP2",
    "EMCN",
    "EPAS1",
    "CDH5"
  ),

  SMC_Pericyte = c(
    "ACTA2",
    "TAGLN",
    "MYH11",
    "CNN1",
    "RGS5",
    "CSPG4",
    "MCAM",
    "PDGFRB"
  ),

  Fibroblast = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "PDGFRA"
  ),

  Myeloid = c(
    "PTPRC",
    "LYZ",
    "LST1",
    "FCER1G",
    "C1QA",
    "C1QB",
    "TYROBP"
  ),

  T_NK = c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "NKG7",
    "GNLY"
  ),

  B_cell = c(
    "MS4A1",
    "CD79A",
    "CD37",
    "CD74"
  ),

  Epithelial = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19"
  )
)

marker_vector <- unique(unlist(marker_sets))
marker_present <- marker_vector[
  marker_vector %in% rownames(obj)
]

cat(
  "\nBroad markers found:",
  length(marker_present),
  "/",
  length(marker_vector),
  "\n"
)

p_dot <- DotPlot(
  obj,
  features = marker_sets,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(PLOT, "10_DotPlot_broad_celltype_markers.png"),
  p_dot,
  width = 18,
  height = 10,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 16. Cluster-level marker summary
# average normalized expression + percent expressing
# ============================================================
cat("\nCalculating marker-panel cluster summary...\n")

rna_data <- LayerData(
  obj,
  assay = "RNA",
  layer = "data"
)

rna_counts <- LayerData(
  obj,
  assay = "RNA",
  layer = "counts"
)

clusters <- as.character(obj$seurat_clusters)

cluster_levels <- sort(
  unique(clusters)
)

marker_summary_list <- list()

idx <- 1

for (cl in cluster_levels) {

  cells_cl <- colnames(obj)[clusters == cl]

  avg_expr <- Matrix::rowMeans(
    rna_data[
      marker_present,
      cells_cl,
      drop = FALSE
    ]
  )

  pct_expr <- Matrix::rowMeans(
    rna_counts[
      marker_present,
      cells_cl,
      drop = FALSE
    ] > 0
  ) * 100

  marker_group <- sapply(
    marker_present,
    function(g) {
      hit <- names(marker_sets)[
        sapply(
          marker_sets,
          function(x) g %in% x
        )
      ]

      if (length(hit) == 0) {
        return("Other")
      }

      paste(hit, collapse = ";")
    }
  )

  marker_summary_list[[idx]] <- data.frame(
    cluster = cl,
    gene = marker_present,
    marker_group = marker_group,
    avg_log_expression = avg_expr,
    percent_expressing = pct_expr,
    n_cells_cluster = length(cells_cl),
    stringsAsFactors = FALSE
  )

  idx <- idx + 1
}

marker_summary <- rbindlist(
  marker_summary_list
)

fwrite(
  marker_summary,
  file.path(
    TAB,
    "08_cluster_marker_panel_summary.csv"
  )
)

group_summary <- marker_summary[
  ,
  .(
    mean_avg_log_expression =
      mean(avg_log_expression),

    mean_percent_expressing =
      mean(percent_expressing)
  ),
  by = .(
    cluster,
    marker_group
  )
]

fwrite(
  group_summary,
  file.path(
    TAB,
    "09_cluster_marker_group_summary.csv"
  )
)


# ============================================================
# 17. Important genes for our study
# Visualization only at Step1
# ============================================================
target_genes <- c(
  "KLF2",
  "KLF4",
  "THBD",
  "EPAS1",
  "ANGPT2",
  "PGF",
  "LOX",
  "CTHRC1",
  "SULF1",
  "SPARC",
  "FN1",
  "COL8A1",
  "MGP",
  "ELN",
  "ITLN1",
  "SLCO4A1"
)

target_present <- target_genes[
  target_genes %in% rownames(obj)
]

writeLines(
  target_present,
  file.path(
    TAB,
    "10_target_genes_present.txt"
  )
)

p_target <- DotPlot(
  obj,
  features = target_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme_bw()

ggsave(
  file.path(PLOT, "11_DotPlot_project_target_genes.png"),
  p_target,
  width = 15,
  height = 8,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 18. Save object
# ============================================================
cat("\nSaving Seurat object...\n")

saveRDS(
  obj,
  file.path(
    OUT,
    "GSE218166_step1_allcells_clustered.rds"
  ),
  compress = FALSE
)


# ============================================================
# 19. Session information
# ============================================================
capture.output(
  sessionInfo(),
  file = file.path(
    LOG,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    paste("Cells:", ncol(obj)),
    paste("Genes:", nrow(obj)),
    paste(
      "Clusters:",
      length(unique(obj$seurat_clusters))
    ),
    paste(
      "FL cells:",
      sum(obj$tech == "FL")
    ),
    paste(
      "3DE cells:",
      sum(obj$tech == "3DE")
    ),
    paste(
      "Normal:",
      sum(obj$disease_stage == "Normal")
    ),
    paste(
      "Early:",
      sum(obj$disease_stage == "Early")
    ),
    paste(
      "Advanced:",
      sum(obj$disease_stage == "Advanced")
    )
  ),
  file.path(
    LOG,
    "step1_summary.txt"
  )
)

cat("\n============================================\n")
cat("STEP1 COMPLETED SUCCESSFULLY\n")
cat("Output:", OUT, "\n")
cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")
cat(
  "Clusters:",
  length(unique(obj$seurat_clusters)),
  "\n"
)
cat("============================================\n")
