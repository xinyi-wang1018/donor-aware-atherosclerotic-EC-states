#!/usr/bin/env Rscript

# ============================================================================
# FINAL strict-source Slingshot refit for atherosclerosis endothelial cells
# Strict sources: GSE155512 + GSE253904 only
# Frozen EC states: C7, C3, C5, C1, C2, C4
# Root: C7
# No forced terminal cluster
# Slingshot input: first 20 dimensions of a locally recomputed Harmony embedding
# ============================================================================

options(stringsAsFactors = FALSE, timeout = 3600)
options(bitmapType = "cairo")
set.seed(20260824)

# ----------------------------- 0. PARAMETERS ---------------------------------
INPUT_RDS <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
DONOR_MAP <- "/root/AS_EC_19day_sprint/04_donor_program_validation/04_sample_to_donor_mapping.csv"
OUTDIR <- "/root/AS_EC_strict_source_slingshot_FINAL"

STRICT_DATASETS <- c("GSE155512", "GSE253904")
CORE_CLUSTERS_NUM <- c("1", "2", "3", "4", "5", "7")
DISPLAY_STATES <- c("C7", "C3", "C5", "C1", "C2", "C4")
MAIN_EXPECTED <- c("C7", "C3", "C5", "C1", "C2")
BRANCH_EXPECTED <- c("C7", "C3", "C4")
ROOT_STATE <- "C7"

EXPECTED_STRICT_CORE_N <- 4435L
N_HVG <- 3000L
N_PCS <- 40L
HARMONY_DIMS <- 1:30
SLINGSHOT_N_DIMS <- 20L
UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.30
N_PT_BINS <- 20L

KEY_GENES <- c(
  "CXCL12", "GJA4", "EFNB2",
  "ANGPT2", "PGF", "LOX", "CTHRC1",
  "SULF1", "ITLN1", "MGP", "ELN", "COL8A1"
)

# -------------------------- 1. PACKAGE CHECK ---------------------------------
required_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "ggplot2", "data.table",
  "harmony", "slingshot", "SingleCellExperiment", "viridisLite"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(data.table)
  library(harmony)
  library(slingshot)
})

# ---------------------------- 2. OUTPUT DIRS ---------------------------------
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
for (d in c("01_objects", "02_plots", "03_tables", "04_logs")) {
  dir.create(file.path(OUTDIR, d), recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(OUTDIR, "04_logs", "strict_source_slingshot_FINAL.log")
log_con <- file(LOG_FILE, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

on.exit({
  cat("\n===== SESSION INFO =====\n")
  print(sessionInfo())
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("===== FINAL STRICT-SOURCE SLINGSHOT REFIT =====\n")
cat("Start time:", format(Sys.time()), "\n")
cat("Input:", INPUT_RDS, "\n")
cat("Strict datasets:", paste(STRICT_DATASETS, collapse = ", "), "\n")
cat("Frozen EC clusters:", paste(CORE_CLUSTERS_NUM, collapse = ", "), "\n")
cat("Root:", ROOT_STATE, "\n")
cat("Forced terminal state: NONE\n")
cat("Output:", OUTDIR, "\n\n")

# ------------------------------ 3. HELPERS -----------------------------------
write_tsv <- function(x, path) {
  data.table::fwrite(as.data.frame(x), path, sep = "\t", na = "NA")
}

save_plot_white <- function(p, path, width = 8, height = 6, dpi = 320) {
  ggplot2::ggsave(
    filename = path, plot = p, width = width, height = height,
    dpi = dpi, bg = "white", limitsize = FALSE
  )
}

standardize_cluster_num <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^[Cc]", "", x)
  x
}

contains_in_order <- function(sequence, required) {
  pos <- match(required, sequence)
  all(!is.na(pos)) && all(diff(pos) > 0)
}

join_rna_layers_if_needed <- function(obj) {
  if (!"RNA" %in% Assays(obj)) stop("RNA assay not found in input object.")
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) {
    obj <- tryCatch(
      SeuratObject::JoinLayers(obj, assay = "RNA"),
      error = function(e) {
        message("JoinLayers skipped: ", conditionMessage(e))
        obj
      }
    )
  }
  obj
}

run_harmony_robust <- function(obj, group_col, dims_use, reduction_save) {
  if (length(unique(as.character(obj@meta.data[[group_col]]))) < 2) {
    stop("Harmony requires >=2 dataset levels after strict filtering.")
  }

  tryCatch(
    harmony::RunHarmony(
      object = obj,
      group.by.vars = group_col,
      reduction.use = "pca",
      dims.use = dims_use,
      reduction.save = reduction_save,
      project.dim = FALSE,
      verbose = TRUE
    ),
    error = function(e_old) {
      message("Legacy Harmony interface failed; retrying current interface: ", conditionMessage(e_old))
      harmony::RunHarmony(
        object = obj,
        group.by.vars = group_col,
        reduction = "pca",
        dims.use = dims_use,
        reduction.save = reduction_save,
        project.dim = FALSE,
        verbose = TRUE
      )
    }
  )
}

fetch_expr <- function(obj, genes) {
  genes_present <- genes[toupper(genes) %in% toupper(rownames(obj))]
  if (length(genes_present) == 0) return(NULL)

  # Preserve actual feature spelling in the object.
  idx <- match(toupper(genes_present), toupper(rownames(obj)))
  genes_present <- unique(rownames(obj)[idx])

  x <- tryCatch(
    Seurat::FetchData(obj, vars = genes_present, layer = "data", clean = FALSE),
    error = function(e) NULL
  )
  if (is.null(x)) {
    x <- Seurat::FetchData(obj, vars = genes_present, slot = "data", clean = FALSE)
  }
  as.data.frame(x)
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(ok)) return(NA_real_)
  if (sum(w[ok]) <= 0) return(mean(x[ok], na.rm = TRUE))
  stats::weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

# ------------------------------- 4. LOAD -------------------------------------
if (!file.exists(INPUT_RDS)) stop("Input RDS not found: ", INPUT_RDS)
obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")

required_meta <- c("dataset_id", "sample_id", "EC_cluster")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Required metadata columns missing: ", paste(missing_meta, collapse = ", "))
}

cat("Loaded input object: ", ncol(obj), " cells; ", nrow(obj), " genes.\n", sep = "")
cat("Datasets in input:\n")
print(table(as.character(obj$dataset_id), useNA = "ifany"))
cat("EC clusters in input:\n")
print(table(as.character(obj$EC_cluster), useNA = "ifany"))

# -------------------------- 5. STRICT FILTER ---------------------------------
obj$dataset_id <- as.character(obj$dataset_id)
obj$sample_id <- as.character(obj$sample_id)
obj$EC_cluster_num <- standardize_cluster_num(obj$EC_cluster)
obj$paper_state <- paste0("C", obj$EC_cluster_num)

keep_cells <- Cells(obj)[
  obj$dataset_id %in% STRICT_DATASETS &
    obj$EC_cluster_num %in% CORE_CLUSTERS_NUM
]

if (length(keep_cells) == 0) stop("Strict filter retained zero cells.")
strict <- subset(obj, cells = keep_cells)
rm(obj); gc()

strict$paper_state <- factor(as.character(strict$paper_state), levels = DISPLAY_STATES)
strict$dataset_id <- factor(as.character(strict$dataset_id), levels = STRICT_DATASETS)

if (!all(STRICT_DATASETS %in% unique(as.character(strict$dataset_id)))) {
  stop("One of the two required strict datasets is absent after filtering.")
}
if (!all(DISPLAY_STATES %in% unique(as.character(strict$paper_state)))) {
  missing_states <- setdiff(DISPLAY_STATES, unique(as.character(strict$paper_state)))
  stop("Required EC state(s) absent after filtering: ", paste(missing_states, collapse = ", "))
}

cat("\nStrict-source core EC retained:", ncol(strict), "\n")
if (ncol(strict) != EXPECTED_STRICT_CORE_N) {
  warning(
    "Strict core cell count is ", ncol(strict),
    ", whereas the current manuscript audit expects ", EXPECTED_STRICT_CORE_N,
    ". Continue, but inspect 01_strict_source_state_counts.tsv before freezing manuscript numbers."
  )
} else {
  cat("Strict core cell count matches expected manuscript audit: ", EXPECTED_STRICT_CORE_N, "\n", sep = "")
}

state_counts <- as.data.frame(
  table(
    dataset_id = as.character(strict$dataset_id),
    paper_state = as.character(strict$paper_state)
  ),
  stringsAsFactors = FALSE
)
names(state_counts)[3] <- "n_cells"
state_counts <- state_counts[state_counts$n_cells > 0, , drop = FALSE]
write_tsv(state_counts, file.path(OUTDIR, "03_tables", "01_strict_source_state_counts.tsv"))

sample_counts <- as.data.frame(
  table(
    dataset_id = as.character(strict$dataset_id),
    sample_id = as.character(strict$sample_id),
    paper_state = as.character(strict$paper_state)
  ),
  stringsAsFactors = FALSE
)
names(sample_counts)[4] <- "n_cells"
sample_counts <- sample_counts[sample_counts$n_cells > 0, , drop = FALSE]
write_tsv(sample_counts, file.path(OUTDIR, "03_tables", "02_strict_source_sample_state_counts.tsv"))

# Optional donor mapping for QC summaries.
strict$donor_id <- NA_character_
if (file.exists(DONOR_MAP)) {
  donor_map <- read.csv(DONOR_MAP, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("dataset_id", "sample_id", "donor_id") %in% colnames(donor_map))) {
    donor_key <- paste(as.character(donor_map$dataset_id), as.character(donor_map$sample_id), sep = "|||MAP|||")
    cell_key <- paste(as.character(strict$dataset_id), as.character(strict$sample_id), sep = "|||MAP|||")
    strict$donor_id <- as.character(donor_map$donor_id[match(cell_key, donor_key)])

    donor_counts <- as.data.frame(
      table(
        dataset_id = as.character(strict$dataset_id),
        donor_id = strict$donor_id,
        paper_state = as.character(strict$paper_state),
        useNA = "ifany"
      ),
      stringsAsFactors = FALSE
    )
    names(donor_counts)[4] <- "n_cells"
    donor_counts <- donor_counts[donor_counts$n_cells > 0, , drop = FALSE]
    write_tsv(donor_counts, file.path(OUTDIR, "03_tables", "03_strict_source_donor_state_counts.tsv"))
  } else {
    warning("Donor map exists but lacks dataset_id/sample_id/donor_id; donor QC table skipped.")
  }
} else {
  message("Donor map not found; donor QC table skipped. This does not block Slingshot.")
}

# -------------------- 6. LOCAL RNA-PCA + HARMONY ----------------------------
cat("\nRecomputing strict-source RNA-PCA/Harmony...\n")
strict <- join_rna_layers_if_needed(strict)
DefaultAssay(strict) <- "RNA"

strict <- NormalizeData(
  strict, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)
strict <- FindVariableFeatures(
  strict, assay = "RNA", selection.method = "vst",
  nfeatures = N_HVG, verbose = FALSE
)
strict <- ScaleData(
  strict, assay = "RNA", features = VariableFeatures(strict), verbose = FALSE
)

max_pcs <- min(N_PCS, length(VariableFeatures(strict)) - 1L, ncol(strict) - 1L)
if (max_pcs < 20) stop("Too few PCs can be computed for the planned analysis: ", max_pcs)

strict <- RunPCA(
  strict, assay = "RNA", features = VariableFeatures(strict),
  npcs = max_pcs, seed.use = 20260824, verbose = FALSE,
  reduction.name = "pca_strict", reduction.key = "PCSTR_"
)

# Harmony wrappers normally expect a reduction called 'pca'.
# Copy the strict PCA to 'pca' so both Harmony interfaces behave consistently.
strict[["pca"]] <- strict[["pca_strict"]]

dims_harmony <- seq_len(min(max(HARMONY_DIMS), max_pcs))
strict <- run_harmony_robust(
  strict, group_col = "dataset_id", dims_use = dims_harmony,
  reduction_save = "harmony_strict"
)

n_harmony_dims <- ncol(Embeddings(strict, "harmony_strict"))
sl_dims <- seq_len(min(SLINGSHOT_N_DIMS, n_harmony_dims))
if (length(sl_dims) < 2) stop("Harmony returned fewer than 2 dimensions.")

strict <- RunUMAP(
  strict,
  reduction = "harmony_strict",
  dims = sl_dims,
  reduction.name = "umap_strict",
  reduction.key = "USTR_",
  n.neighbors = UMAP_N_NEIGHBORS,
  min.dist = UMAP_MIN_DIST,
  metric = "cosine",
  seed.use = 20260824,
  verbose = FALSE
)

reduction_info <- data.frame(
  item = c(
    "n_cells", "n_genes", "n_variable_features", "n_pca",
    "n_harmony", "harmony_dims_for_slingshot", "root", "forced_end"
  ),
  value = c(
    ncol(strict), nrow(strict), length(VariableFeatures(strict)), max_pcs,
    n_harmony_dims, paste(sl_dims, collapse = ","), ROOT_STATE, "NONE"
  ),
  stringsAsFactors = FALSE
)
write_tsv(reduction_info, file.path(OUTDIR, "03_tables", "04_reduction_and_model_parameters.tsv"))

# ---------------------------- 7. SLINGSHOT -----------------------------------
cat("\nRunning Slingshot on strict Harmony dimensions: ", paste(sl_dims, collapse = ","), "\n", sep = "")
emb <- Embeddings(strict, "harmony_strict")[, sl_dims, drop = FALSE]
labels <- droplevels(factor(as.character(strict$paper_state), levels = DISPLAY_STATES))

if (!ROOT_STATE %in% levels(labels)) stop("Root state ", ROOT_STATE, " is absent.")

sds <- slingshot::slingshot(
  data = as.matrix(emb),
  clusterLabels = labels,
  start.clus = ROOT_STATE,
  stretch = 0,
  reweight = TRUE,
  reassign = TRUE
)

lineages <- slingshot::slingLineages(sds)
if (length(lineages) == 0) stop("Slingshot returned no lineage.")

lineage_names <- names(lineages)
if (is.null(lineage_names) || any(lineage_names == "")) {
  lineage_names <- paste0("Lineage", seq_along(lineages))
}

lineage_table <- do.call(rbind, lapply(seq_along(lineages), function(i) {
  seq_i <- as.character(lineages[[i]])
  n_main <- sum(MAIN_EXPECTED %in% seq_i)
  score <- 10 * n_main - 5 * as.integer("C4" %in% seq_i)
  if (contains_in_order(seq_i, MAIN_EXPECTED)) score <- score + 100
  data.frame(
    lineage_index = i,
    lineage = lineage_names[i],
    sequence = paste(seq_i, collapse = "->"),
    endpoint = tail(seq_i, 1),
    exact_main_C7_C3_C5_C1_C2 = contains_in_order(seq_i, MAIN_EXPECTED),
    branch_C7_C3_C4 = contains_in_order(seq_i, BRANCH_EXPECTED),
    n_main_states_present = n_main,
    primary_score = score,
    stringsAsFactors = FALSE
  )
}))

primary_idx <- which.max(lineage_table$primary_score)
lineage_table$is_primary <- seq_len(nrow(lineage_table)) == primary_idx
write_tsv(lineage_table, file.path(OUTDIR, "03_tables", "05_slingshot_lineages.tsv"))

cat("Slingshot lineage sequences:\n")
print(lineage_table)
cat("Selected primary lineage:", lineage_table$lineage[primary_idx], "=", lineage_table$sequence[primary_idx], "\n")

pt <- slingshot::slingPseudotime(sds, na = TRUE)
wt <- slingshot::slingCurveWeights(sds)
if (is.null(dim(pt))) pt <- matrix(pt, ncol = 1)
if (is.null(dim(wt))) wt <- matrix(wt, ncol = 1)
rownames(pt) <- rownames(emb)
rownames(wt) <- rownames(emb)

if (ncol(pt) != length(lineages)) {
  warning("Pseudotime column count differs from lineage count; verify lineage tables manually.")
}

pt_colnames <- paste0("strict_pseudotime_", lineage_names[seq_len(ncol(pt))])
wt_colnames <- paste0("strict_weight_", lineage_names[seq_len(ncol(wt))])
colnames(pt) <- pt_colnames
colnames(wt) <- wt_colnames

for (j in seq_len(ncol(pt))) strict[[pt_colnames[j]]] <- pt[Cells(strict), j]
for (j in seq_len(ncol(wt))) strict[[wt_colnames[j]]] <- wt[Cells(strict), j]
strict$strict_primary_lineage <- lineage_names[primary_idx]
strict$strict_primary_pseudotime <- pt[Cells(strict), primary_idx]
strict$strict_primary_weight <- wt[Cells(strict), primary_idx]

# Save Slingshot object now.
saveRDS(sds, file.path(OUTDIR, "01_objects", "strict_source_slingshot_sds_FINAL.rds"), compress = TRUE)

# ---------------------- 8. CELL-LEVEL OUTPUT ---------------------------------
um <- Embeddings(strict, "umap_strict")
cell_out <- data.frame(
  cell_id = Cells(strict),
  UMAP_1 = um[Cells(strict), 1],
  UMAP_2 = um[Cells(strict), 2],
  dataset_id = as.character(strict$dataset_id),
  sample_id = as.character(strict$sample_id),
  donor_id = strict$donor_id,
  EC_cluster = as.character(strict$EC_cluster_num),
  paper_state = as.character(strict$paper_state),
  primary_lineage = strict$strict_primary_lineage,
  primary_pseudotime = strict$strict_primary_pseudotime,
  primary_weight = strict$strict_primary_weight,
  stringsAsFactors = FALSE
)

for (j in seq_len(ncol(pt))) cell_out[[pt_colnames[j]]] <- pt[Cells(strict), j]
for (j in seq_len(ncol(wt))) cell_out[[wt_colnames[j]]] <- wt[Cells(strict), j]
write_tsv(cell_out, file.path(OUTDIR, "03_tables", "06_cell_level_strict_source_pseudotime_weights.tsv"))

# ------------------- 9. STATE PSEUDOTIME SUMMARY ----------------------------
state_summary_all <- do.call(rbind, lapply(seq_len(ncol(pt)), function(j) {
  do.call(rbind, lapply(DISPLAY_STATES, function(st) {
    ids <- which(as.character(strict$paper_state) == st)
    v <- pt[ids, j]
    w <- wt[ids, j]
    vf <- v[is.finite(v)]
    wf <- w[is.finite(w)]
    data.frame(
      lineage = lineage_names[j],
      paper_state = st,
      n_cells_total = length(ids),
      n_cells_finite_pseudotime = length(vf),
      median_pseudotime = if (length(vf)) median(vf) else NA_real_,
      q25_pseudotime = if (length(vf)) unname(quantile(vf, 0.25)) else NA_real_,
      q75_pseudotime = if (length(vf)) unname(quantile(vf, 0.75)) else NA_real_,
      median_lineage_weight = if (length(wf)) median(wf) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))
write_tsv(state_summary_all, file.path(OUTDIR, "03_tables", "07_state_pseudotime_summary_all_lineages.tsv"))

primary_name <- lineage_names[primary_idx]
primary_summary <- state_summary_all[state_summary_all$lineage == primary_name, , drop = FALSE]
write_tsv(primary_summary, file.path(OUTDIR, "03_tables", "08_PRIMARY_state_pseudotime_summary.tsv"))

main_primary <- primary_summary[primary_summary$paper_state %in% MAIN_EXPECTED, , drop = FALSE]
main_primary$expected_rank <- match(main_primary$paper_state, MAIN_EXPECTED)
main_primary$observed_rank <- rank(main_primary$median_pseudotime, ties.method = "min", na.last = "keep")
main_primary <- main_primary[order(main_primary$expected_rank), , drop = FALSE]
main_primary$order_agrees <- main_primary$expected_rank == main_primary$observed_rank
write_tsv(main_primary, file.path(OUTDIR, "03_tables", "09_PRIMARY_expected_vs_observed_state_order.tsv"))

observed_main_order <- main_primary$paper_state[order(main_primary$median_pseudotime, na.last = NA)]
exact_median_order <- identical(observed_main_order, MAIN_EXPECTED)

# ---------------------- 10. GENE EXPRESSION ----------------------------------
expr <- fetch_expr(strict, KEY_GENES)
if (!is.null(expr)) {
  genes_present <- colnames(expr)

  # Cell-level key-gene output linked to primary pseudotime.
  gene_cell <- data.frame(
    cell_id = rownames(expr),
    paper_state = as.character(strict$paper_state[match(rownames(expr), Cells(strict))]),
    primary_pseudotime = strict$strict_primary_pseudotime[match(rownames(expr), Cells(strict))],
    primary_weight = strict$strict_primary_weight[match(rownames(expr), Cells(strict))],
    expr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write_tsv(gene_cell, file.path(OUTDIR, "03_tables", "10_key_gene_cell_expression_primary_pseudotime.tsv"))

  # State-level expression summary.
  state_gene_list <- list()
  k <- 1L
  for (st in DISPLAY_STATES) {
    ids <- which(gene_cell$paper_state == st)
    if (!length(ids)) next
    for (g in genes_present) {
      v <- as.numeric(gene_cell[ids, g])
      state_gene_list[[k]] <- data.frame(
        paper_state = st,
        gene = g,
        n_cells = length(v),
        mean_logexpr = mean(v, na.rm = TRUE),
        median_logexpr = median(v, na.rm = TRUE),
        pct_expressing = mean(v > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  state_gene_summary <- do.call(rbind, state_gene_list)
  write_tsv(state_gene_summary, file.path(OUTDIR, "03_tables", "11_key_gene_expression_by_state.tsv"))

  # Descriptive Spearman correlation with primary pseudotime (NO inferential p-value).
  corr_list <- lapply(genes_present, function(g) {
    ok <- is.finite(gene_cell$primary_pseudotime) & is.finite(gene_cell$primary_weight)
    if (sum(ok) < 20) {
      rho <- NA_real_
    } else {
      rho <- suppressWarnings(cor(
        gene_cell$primary_pseudotime[ok], as.numeric(gene_cell[ok, g]),
        method = "spearman", use = "complete.obs"
      ))
    }
    data.frame(gene = g, n_cells = sum(ok), spearman_rho_descriptive = rho, stringsAsFactors = FALSE)
  })
  corr_df <- do.call(rbind, corr_list)
  write_tsv(corr_df, file.path(OUTDIR, "03_tables", "12_key_gene_pseudotime_descriptive_correlations.tsv"))

  # Pseudotime bins with lineage-weighted mean expression.
  ok_pt <- is.finite(gene_cell$primary_pseudotime)
  ranks <- rank(gene_cell$primary_pseudotime[ok_pt], ties.method = "average")
  bins <- pmin(N_PT_BINS, pmax(1L, ceiling(ranks / max(ranks) * N_PT_BINS)))
  gene_cell$pt_bin <- NA_integer_
  gene_cell$pt_bin[ok_pt] <- bins

  bin_list <- list()
  z <- 1L
  for (b in sort(unique(na.omit(gene_cell$pt_bin)))) {
    ids <- which(gene_cell$pt_bin == b)
    for (g in genes_present) {
      bin_list[[z]] <- data.frame(
        pt_bin = b,
        n_cells = length(ids),
        median_pseudotime = median(gene_cell$primary_pseudotime[ids], na.rm = TRUE),
        gene = g,
        weighted_mean_logexpr = weighted_mean_safe(
          as.numeric(gene_cell[ids, g]),
          pmax(gene_cell$primary_weight[ids], 0)
        ),
        stringsAsFactors = FALSE
      )
      z <- z + 1L
    }
  }
  bin_df <- do.call(rbind, bin_list)
  write_tsv(bin_df, file.path(OUTDIR, "03_tables", "13_key_gene_primary_pseudotime_binned_trends.tsv"))
}

# ------------------------------- 11. PLOTS -----------------------------------
white_theme <- theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = "black")
  )

state_df <- data.frame(
  UMAP_1 = um[, 1], UMAP_2 = um[, 2],
  paper_state = factor(as.character(strict$paper_state), levels = DISPLAY_STATES),
  dataset_id = as.character(strict$dataset_id),
  primary_pseudotime = strict$strict_primary_pseudotime,
  primary_weight = strict$strict_primary_weight,
  stringsAsFactors = FALSE
)

p_state <- ggplot(state_df, aes(UMAP_1, UMAP_2, color = paper_state)) +
  geom_point(size = 0.75, alpha = 0.80) +
  white_theme +
  labs(title = "Strict-source core EC states", color = "EC state")
save_plot_white(p_state, file.path(OUTDIR, "02_plots", "01_strict_source_state_UMAP.png"), 7.5, 6.2)

p_dataset <- ggplot(state_df, aes(UMAP_1, UMAP_2, color = dataset_id)) +
  geom_point(size = 0.75, alpha = 0.75) +
  white_theme +
  labs(title = "Strict-source dataset distribution", color = "Dataset")
save_plot_white(p_dataset, file.path(OUTDIR, "02_plots", "02_strict_source_dataset_UMAP.png"), 7.5, 6.2)

p_pt <- ggplot(state_df, aes(UMAP_1, UMAP_2, color = primary_pseudotime)) +
  geom_point(size = 0.75, alpha = 0.82) +
  scale_color_gradientn(colors = viridisLite::viridis(256), na.value = "grey88") +
  white_theme +
  labs(
    title = paste0("Primary Slingshot pseudotime: ", lineage_table$sequence[primary_idx]),
    color = "Pseudotime"
  )
save_plot_white(p_pt, file.path(OUTDIR, "02_plots", "03_PRIMARY_pseudotime_UMAP.png"), 8.4, 6.5)

p_wt <- ggplot(state_df, aes(UMAP_1, UMAP_2, color = primary_weight)) +
  geom_point(size = 0.75, alpha = 0.82) +
  scale_color_gradientn(colors = viridisLite::viridis(256), na.value = "grey88") +
  white_theme +
  labs(title = "Primary-lineage curve weight", color = "Weight")
save_plot_white(p_wt, file.path(OUTDIR, "02_plots", "04_PRIMARY_lineage_weight_UMAP.png"), 7.5, 6.2)

box_df <- data.frame(
  paper_state = factor(as.character(strict$paper_state), levels = DISPLAY_STATES),
  primary_pseudotime = strict$strict_primary_pseudotime
)
p_box <- ggplot(box_df[is.finite(box_df$primary_pseudotime), ], aes(paper_state, primary_pseudotime)) +
  geom_boxplot(outlier.size = 0.25, width = 0.65) +
  white_theme +
  labs(title = "Primary-lineage pseudotime by frozen EC state", x = "EC state", y = "Pseudotime")
save_plot_white(p_box, file.path(OUTDIR, "02_plots", "05_PRIMARY_state_pseudotime_boxplot.png"), 8.2, 5.8)

# Per-lineage pseudotime UMAPs.
for (j in seq_len(ncol(pt))) {
  d <- state_df
  d$pt_lineage <- pt[Cells(strict), j]
  pp <- ggplot(d, aes(UMAP_1, UMAP_2, color = pt_lineage)) +
    geom_point(size = 0.70, alpha = 0.82) +
    scale_color_gradientn(colors = viridisLite::viridis(256), na.value = "grey88") +
    white_theme +
    labs(
      title = paste0(lineage_names[j], ": ", lineage_table$sequence[j]),
      color = "Pseudotime"
    )
  save_plot_white(
    pp,
    file.path(OUTDIR, "02_plots", sprintf("06_lineage_%02d_pseudotime_UMAP.png", j)),
    8.2, 6.2
  )
}

# Actual Slingshot principal curves in the first two strict Harmony dimensions.
# This plot is a geometry/QC figure; pseudotime inference itself used all sl_dims above.
curve_png <- file.path(OUTDIR, "02_plots", "07_slingshot_curves_Harmony12.png")
png(curve_png, width = 2400, height = 1900, res = 280, bg = "white")
tryCatch({
  st_chr <- as.character(strict$paper_state)
  st_cols <- setNames(grDevices::hcl.colors(length(DISPLAY_STATES), "Dark 3"), DISPLAY_STATES)
  plot(
    emb[, 1], emb[, 2],
    col = grDevices::adjustcolor(st_cols[st_chr], alpha.f = 0.62),
    pch = 16, cex = 0.55,
    xlab = "Harmony 1", ylab = "Harmony 2",
    main = "Strict-source Slingshot principal curves"
  )
  try(lines(sds, lwd = 4, col = "black"), silent = TRUE)
  legend(
    "topright", legend = DISPLAY_STATES, col = st_cols[DISPLAY_STATES],
    pch = 16, bty = "n", title = "EC state", cex = 0.85
  )
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, paste("Curve plotting failed but Slingshot model succeeded:\n", conditionMessage(e)))
})
dev.off()

# Binned gene-trend plot if available.
if (exists("bin_df") && nrow(bin_df) > 0) {
  p_gene <- ggplot(
    bin_df,
    aes(x = median_pseudotime, y = weighted_mean_logexpr, group = gene)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    white_theme +
    labs(
      title = "Key-gene trends along the primary strict-source pseudotime",
      x = "Median pseudotime bin",
      y = "Lineage-weighted mean log-normalized expression"
    )
  save_plot_white(p_gene, file.path(OUTDIR, "02_plots", "08_PRIMARY_key_gene_pseudotime_trends.png"), 12.5, 8.0)
}

# -------------------------- 12. FINAL SUMMARY --------------------------------
primary_sequence <- lineage_table$sequence[primary_idx]
main_exact_detected <- lineage_table$exact_main_C7_C3_C5_C1_C2[primary_idx]
branch_detected <- any(lineage_table$branch_C7_C3_C4)

summary_lines <- c(
  "FINAL STRICT-SOURCE SLINGSHOT REFIT",
  paste0("Input: ", INPUT_RDS),
  paste0("Strict datasets: ", paste(STRICT_DATASETS, collapse = " + ")),
  paste0("Strict core EC cells: ", ncol(strict)),
  paste0("Expected strict core EC cells from manuscript audit: ", EXPECTED_STRICT_CORE_N),
  paste0("Primary lineage: ", primary_sequence),
  paste0("Exact C7->C3->C5->C1->C2 sequence recovered: ", main_exact_detected),
  paste0("C7->C3->C4 branch recovered in any lineage: ", branch_detected),
  paste0("Observed primary-lineage state median order: ", paste(observed_main_order, collapse = "->")),
  paste0("Median pseudotime order exactly matches C7->C3->C5->C1->C2: ", exact_median_order),
  "",
  "IMPORTANT INTERPRETATION:",
  "Pseudotime is a relative transcriptional organization, not proof of real lineage conversion or chronological disease progression.",
  "No terminal state was forced; only C7 was specified as the root.",
  ""
)

writeLines(summary_lines, file.path(OUTDIR, "FINAL_RESULT_SUMMARY.txt"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")

# Save final strict object with all new metadata/reductions.
saveRDS(
  strict,
  file.path(OUTDIR, "01_objects", "AS_EC_strict_source_core_slingshot_FINAL.rds"),
  compress = TRUE
)

# Lightweight manifest.
manifest <- data.frame(
  file = list.files(OUTDIR, recursive = TRUE, full.names = FALSE),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(OUTDIR, "03_tables", "99_output_manifest.tsv"))

cat("\n===== SUCCESS =====\n")
cat("End time:", format(Sys.time()), "\n")
cat("Primary result summary:", file.path(OUTDIR, "FINAL_RESULT_SUMMARY.txt"), "\n")
cat("Final object:", file.path(OUTDIR, "01_objects", "AS_EC_strict_source_core_slingshot_FINAL.rds"), "\n")
cat("Upload/review package should include FINAL_RESULT_SUMMARY.txt + 02_plots + 03_tables + 04_logs.\n")
