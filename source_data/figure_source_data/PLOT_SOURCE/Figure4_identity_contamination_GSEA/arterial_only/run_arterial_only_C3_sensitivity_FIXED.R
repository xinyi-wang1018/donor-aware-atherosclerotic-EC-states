#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260812)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

# =========================
# 0. Paths and parameters
# =========================
input_rds <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
outdir <- "/root/AS_EC_arterial_only_C3_sensitivity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "01_objects"), showWarnings = FALSE)
dir.create(file.path(outdir, "02_plots"), showWarnings = FALSE)
dir.create(file.path(outdir, "03_tables"), showWarnings = FALSE)
dir.create(file.path(outdir, "04_logs"), showWarnings = FALSE)

log_file <- file.path(outdir, "04_logs", "arterial_only_C3_sensitivity.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("===== STRICT ARTERIAL-ONLY C3 SENSITIVITY ANALYSIS =====\n")
cat("Start time:", format(Sys.time()), "\n")
cat("Input:", input_rds, "\n")
cat("Output:", outdir, "\n\n")

if (!file.exists(input_rds)) {
  stop("Input RDS not found: ", input_rds)
}
if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("Package 'harmony' is required. It should already exist because the project previously used Harmony.")
}
if (!requireNamespace("slingshot", quietly = TRUE)) {
  stop("Package 'slingshot' is required. It should already exist because the project previously used Slingshot.")
}

# =========================
# 1. Helper functions
# =========================
white_theme <- theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = "black")
  )

save_plot_white <- function(p, filename, width = 8, height = 6, dpi = 300) {
  ggsave(filename = filename, plot = p, width = width, height = height,
         dpi = dpi, bg = "white")
}

canonical_state <- function(x) {
  # Accept both explicit labels (C0-C7) and numeric EC_cluster labels (0-7).
  # This object stores the original EC cluster mainly in EC_cluster as numeric values.
  raw <- trimws(as.character(x))
  up <- toupper(raw)
  out <- rep(NA_character_, length(up))

  # 1) Explicit C0-C7 anywhere in the label.
  has_c <- grepl("C\\s*[0-7]([^0-9]|$)", up, perl = TRUE)
  if (any(has_c)) {
    tmp <- sub("^.*?C\\s*([0-7])([^0-9]|$).*$", "C\\1", up[has_c], perl = TRUE)
    out[has_c] <- tmp
  }

  # 2) Bare numeric cluster labels 0-7, e.g. EC_cluster = 3 -> C3.
  bare_num <- grepl("^[0-7]$", raw) & is.na(out)
  out[bare_num] <- paste0("C", raw[bare_num])

  out
}

detect_state_col <- function(md) {
  preferred <- c(
    "EC_state", "ec_state", "State", "state", "EC_cluster", "ec_cluster",
    "EC_subcluster", "ec_subcluster", "cell_state", "CellState",
    "subcluster", "sub_cluster", "cluster_label", "Cluster"
  )
  candidates <- c(intersect(preferred, colnames(md)), setdiff(colnames(md), preferred))
  best_col <- NULL
  best_score <- -Inf
  for (nm in candidates) {
    v <- canonical_state(md[[nm]])
    frac <- mean(!is.na(v))
    nstates <- length(unique(na.omit(v)))
    pref_bonus <- ifelse(nm %in% preferred, 2, 0)
    score <- 10 * frac + nstates + pref_bonus
    if (nstates >= 3 && score > best_score) {
      best_score <- score
      best_col <- nm
    }
  }
  if (is.null(best_col)) {
    stop("Could not auto-detect a metadata column containing C0-C7 state labels. Metadata columns are: ",
         paste(colnames(md), collapse = ", "))
  }
  best_col
}

detect_meta_col <- function(md, preferred, min_unique = 2) {
  for (nm in preferred) {
    if (nm %in% colnames(md)) {
      u <- unique(na.omit(as.character(md[[nm]])))
      if (length(u) >= min_unique) return(nm)
    }
  }
  NULL
}

match_genes_ci <- function(wanted, features) {
  idx <- match(toupper(wanted), toupper(features))
  out <- features[idx[!is.na(idx)]]
  unique(out)
}

fetch_data_layer <- function(obj, vars) {
  if (length(vars) == 0) return(data.frame(row.names = colnames(obj)))
  x <- tryCatch(
    FetchData(obj, vars = vars, layer = "data", clean = FALSE),
    error = function(e) NULL
  )
  if (is.null(x)) {
    x <- FetchData(obj, vars = vars, slot = "data", clean = FALSE)
  }
  as.data.frame(x)
}

program_zscore <- function(obj, genes) {
  genes2 <- match_genes_ci(genes, rownames(obj))
  if (length(genes2) == 0) {
    return(list(score = rep(NA_real_, ncol(obj)), genes = character(0)))
  }
  x <- fetch_data_layer(obj, genes2)
  zmat <- sapply(x, function(v) {
    v <- as.numeric(v)
    s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / s
  })
  if (is.null(dim(zmat))) zmat <- matrix(zmat, ncol = 1)
  score <- rowMeans(zmat, na.rm = TRUE)
  names(score) <- rownames(x)
  list(score = score, genes = genes2)
}

safe_quantile <- function(x, p = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 10) stop("Too few finite values to estimate a strict threshold.")
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# =========================
# 2. Load object and identify metadata
# =========================
obj <- readRDS(input_rds)
cat("Loaded object with", ncol(obj), "cells and", nrow(obj), "features.\n")

if (!("RNA" %in% Assays(obj))) stop("RNA assay not found in input object.")
DefaultAssay(obj) <- "RNA"

# Seurat v5: join split RNA layers if needed, making downstream DE/fetch robust.
if (inherits(obj[["RNA"]], "Assay5") && "JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj <- tryCatch(
    SeuratObject::JoinLayers(obj, assay = "RNA"),
    error = function(e) {
      cat("JoinLayers skipped because:", conditionMessage(e), "\n")
      obj
    }
  )
}

# Ensure log-normalized RNA data are available for signature scoring and DE.
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

md <- obj@meta.data
state_col <- detect_state_col(md)
dataset_col <- detect_meta_col(md, c("dataset_id", "dataset", "Dataset", "cohort", "cohort_id", "study", "Study", "GEO", "geo"))
sample_col <- detect_meta_col(md, c("sample_id", "sample", "Sample", "sample_name", "orig.ident", "donor_id", "donor", "patient_id", "patient", "subject_id", "subject"))

obj$arterial_state_original <- canonical_state(md[[state_col]])

if (is.null(sample_col)) {
  obj$detected_sample <- "all_cells"
  sample_col_use <- "detected_sample"
} else {
  obj$detected_sample <- as.character(obj@meta.data[[sample_col]])
  sample_col_use <- "detected_sample"
}
if (is.null(dataset_col)) {
  # If no explicit dataset column exists, use sample as the Harmony grouping variable.
  obj$detected_dataset <- obj$detected_sample
  dataset_col_use <- "detected_dataset"
} else {
  obj$detected_dataset <- as.character(obj@meta.data[[dataset_col]])
  dataset_col_use <- "detected_dataset"
}

meta_detect <- data.frame(
  role = c("state", "dataset", "sample"),
  detected_column = c(state_col,
                      ifelse(is.null(dataset_col), "[not found; using sample]", dataset_col),
                      ifelse(is.null(sample_col), "[not found; all_cells]", sample_col)),
  stringsAsFactors = FALSE
)
write_tsv(meta_detect, file.path(outdir, "03_tables", "00_detected_metadata_columns.tsv"))
cat("Detected state column:", state_col, "\n")
cat("Detected dataset column:", ifelse(is.null(dataset_col), "NONE", dataset_col), "\n")
cat("Detected sample column:", ifelse(is.null(sample_col), "NONE", sample_col), "\n\n")

# =========================
# 3. Marker programs
# =========================
# Core vascular-type markers used for contamination exclusion.
marker_sets <- list(
  arterial_core = c("GJA4", "GJA5", "SEMA3G"),
  arterial_project = c("GJA4", "GJA5", "SEMA3G", "EFNB2", "SOX17", "HEY1", "CXCL12"),
  venous_core = c("ACKR1", "NR2F2"),
  capillary_core = c("RGCC", "CA4")
)

marker_availability <- do.call(rbind, lapply(names(marker_sets), function(nm) {
  g <- marker_sets[[nm]]
  present <- match_genes_ci(g, rownames(obj))
  data.frame(program = nm, requested = paste(g, collapse = ","),
             detected = paste(present, collapse = ","),
             n_requested = length(g), n_detected = length(present),
             stringsAsFactors = FALSE)
}))
write_tsv(marker_availability, file.path(outdir, "03_tables", "01_marker_sets_used.tsv"))

if (marker_availability$n_detected[marker_availability$program == "venous_core"] < 1 ||
    marker_availability$n_detected[marker_availability$program == "capillary_core"] < 1 ||
    marker_availability$n_detected[marker_availability$program == "arterial_project"] < 2) {
  stop("Insufficient vascular-type marker genes detected in the RNA assay. See 01_marker_sets_used.tsv")
}

sc_art <- program_zscore(obj, marker_sets$arterial_project)
sc_ven <- program_zscore(obj, marker_sets$venous_core)
sc_cap <- program_zscore(obj, marker_sets$capillary_core)

# Match scores back by cell name.
obj$arterial_program_z <- sc_art$score[colnames(obj)]
obj$venous_program_z <- sc_ven$score[colnames(obj)]
obj$capillary_program_z <- sc_cap$score[colnames(obj)]
obj$nonarterial_program_z <- pmax(obj$venous_program_z, obj$capillary_program_z, na.rm = TRUE)

axis_states <- c("C7", "C3", "C5", "C1", "C2")
state_factor_levels <- c("C7", "C3", "C5", "C1", "C2")

# =========================
# 4. Define strict arterial-only cells
# =========================
# Use C7 (known arterial-stable root) as the negative-control background for
# venous/capillary programs. A cell is removed if its venous or capillary score
# exceeds the 95th percentile seen in C7.
idx_c7 <- obj$arterial_state_original == "C7"
if (sum(idx_c7, na.rm = TRUE) < 50) stop("Fewer than 50 C7 cells found; cannot define robust C7-based thresholds.")

ven_cut <- safe_quantile(obj$venous_program_z[idx_c7], 0.95)
cap_cut <- safe_quantile(obj$capillary_program_z[idx_c7], 0.95)

obj$axis_state_flag <- obj$arterial_state_original %in% axis_states
obj$primary_arterial_only <- obj$axis_state_flag &
  is.finite(obj$venous_program_z) & is.finite(obj$capillary_program_z) &
  obj$venous_program_z <= ven_cut & obj$capillary_program_z <= cap_cut

# Secondary, deliberately more stringent sensitivity layer:
# require the arterial program to dominate both non-arterial programs.
obj$ultra_arterial_dominant <- obj$primary_arterial_only &
  is.finite(obj$arterial_program_z) &
  obj$arterial_program_z > obj$venous_program_z &
  obj$arterial_program_z > obj$capillary_program_z

thresholds <- data.frame(
  threshold = c("C7_venous_program_95pct", "C7_capillary_program_95pct"),
  value = c(ven_cut, cap_cut)
)
write_tsv(thresholds, file.path(outdir, "03_tables", "02_filter_thresholds.tsv"))

count_state <- function(flag_name) {
  st <- obj$arterial_state_original
  before <- as.data.frame(table(factor(st[st %in% axis_states], levels = state_factor_levels)))
  colnames(before) <- c("state", "n_before")
  keep <- obj@meta.data[[flag_name]] %in% TRUE
  after <- as.data.frame(table(factor(st[keep & st %in% axis_states], levels = state_factor_levels)))
  colnames(after) <- c("state", "n_after")
  out <- merge(before, after, by = "state", all = TRUE)
  out$retention_fraction <- ifelse(out$n_before > 0, out$n_after / out$n_before, NA_real_)
  out$filter <- flag_name
  out[, c("filter", "state", "n_before", "n_after", "retention_fraction")]
}

state_counts <- rbind(count_state("primary_arterial_only"), count_state("ultra_arterial_dominant"))
write_tsv(state_counts, file.path(outdir, "03_tables", "03_state_cell_counts_before_after.tsv"))

# C3 coverage across datasets and samples after primary filter.
md2 <- obj@meta.data
c3_cov <- md2[md2$arterial_state_original == "C3", c("detected_dataset", "detected_sample", "primary_arterial_only", "ultra_arterial_dominant"), drop = FALSE]
if (nrow(c3_cov) > 0) {
  c3_cov$n <- 1L
  cov_ag <- aggregate(n ~ detected_dataset + detected_sample + primary_arterial_only + ultra_arterial_dominant,
                      data = c3_cov, FUN = sum)
  write_tsv(cov_ag, file.path(outdir, "03_tables", "04_C3_dataset_sample_coverage.tsv"))
}

# Score summaries by original state.
score_df <- obj@meta.data[obj$arterial_state_original %in% axis_states,
                          c("arterial_state_original", "arterial_program_z", "venous_program_z", "capillary_program_z"), drop = FALSE]
score_summary <- do.call(rbind, lapply(split(score_df, score_df$arterial_state_original), function(d) {
  data.frame(
    state = d$arterial_state_original[1],
    n = nrow(d),
    arterial_median = median(d$arterial_program_z, na.rm = TRUE),
    venous_median = median(d$venous_program_z, na.rm = TRUE),
    capillary_median = median(d$capillary_program_z, na.rm = TRUE),
    arterial_q25 = quantile(d$arterial_program_z, 0.25, na.rm = TRUE),
    arterial_q75 = quantile(d$arterial_program_z, 0.75, na.rm = TRUE),
    venous_q95 = quantile(d$venous_program_z, 0.95, na.rm = TRUE),
    capillary_q95 = quantile(d$capillary_program_z, 0.95, na.rm = TRUE)
  )
}))
score_summary$state <- factor(score_summary$state, levels = state_factor_levels)
score_summary <- score_summary[order(score_summary$state), ]
write_tsv(score_summary, file.path(outdir, "03_tables", "05_program_scores_by_state.tsv"))

# Score violin plot before filtering.
score_long <- rbind(
  data.frame(state = score_df$arterial_state_original, program = "Arterial", score = score_df$arterial_program_z),
  data.frame(state = score_df$arterial_state_original, program = "Venous", score = score_df$venous_program_z),
  data.frame(state = score_df$arterial_state_original, program = "Capillary", score = score_df$capillary_program_z)
)
score_long$state <- factor(score_long$state, levels = state_factor_levels)
p_scores <- ggplot(score_long, aes(x = state, y = score)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  facet_wrap(~ program, scales = "free_y", nrow = 1) +
  labs(x = "Original EC state", y = "Standardized program score",
       title = "Vascular identity programs before strict arterial-only filtering") +
  white_theme
save_plot_white(p_scores, file.path(outdir, "02_plots", "Fig01_vascular_program_scores_before_filter.png"), 12, 4.5)

# =========================
# 5. Subset primary strict arterial-only cells
# =========================
keep_cells <- colnames(obj)[obj$primary_arterial_only %in% TRUE]
if (length(keep_cells) < 200) stop("Primary arterial-only filter retained fewer than 200 cells.")
obj_strict <- subset(obj, cells = keep_cells)

# Preserve the scores and labels; redo feature selection/PCA/Harmony/UMAP/clustering from the restricted cells.
DefaultAssay(obj_strict) <- "RNA"
obj_strict <- NormalizeData(obj_strict, verbose = FALSE)
obj_strict <- FindVariableFeatures(obj_strict, selection.method = "vst",
                                   nfeatures = min(3000, max(500, nrow(obj_strict) - 1)),
                                   verbose = FALSE)
obj_strict <- ScaleData(obj_strict, features = VariableFeatures(obj_strict), verbose = FALSE)

npcs <- min(30, length(VariableFeatures(obj_strict)) - 1, ncol(obj_strict) - 1)
if (npcs < 10) stop("Too few dimensions available after filtering.")
obj_strict <- RunPCA(obj_strict, features = VariableFeatures(obj_strict), npcs = npcs,
                     verbose = FALSE)

harmony_group <- dataset_col_use
n_groups <- length(unique(na.omit(as.character(obj_strict@meta.data[[harmony_group]]))))
if (n_groups >= 2) {
  obj_strict <- harmony::RunHarmony(
    obj_strict,
    group.by.vars = harmony_group,
    reduction = "pca",
    dims.use = seq_len(npcs),
    assay.use = "RNA",
    reduction.save = "harmony_arterial_only",
    verbose = FALSE
  )
  reduction_use <- "harmony_arterial_only"
} else {
  cat("Only one Harmony group detected; using PCA without Harmony.\n")
  reduction_use <- "pca"
}

dims_use <- seq_len(min(20, ncol(Embeddings(obj_strict, reduction = reduction_use))))
obj_strict <- FindNeighbors(obj_strict, reduction = reduction_use, dims = dims_use, verbose = FALSE)
obj_strict <- FindClusters(obj_strict, resolution = 0.4, random.seed = 20260812, verbose = FALSE)
obj_strict$strict_recluster <- as.character(Idents(obj_strict))
obj_strict <- RunUMAP(obj_strict, reduction = reduction_use, dims = dims_use,
                      reduction.name = "umap_arterial_only", seed.use = 20260812,
                      verbose = FALSE)

# Save strict object before pseudotime.
saveRDS(obj_strict, file.path(outdir, "01_objects", "AS_EC_strict_arterial_only_reclustered.rds"), compress = FALSE)

# UMAP by original state and by new strict cluster.
p_state <- DimPlot(obj_strict, reduction = "umap_arterial_only",
                   group.by = "arterial_state_original", label = TRUE, repel = TRUE) +
  labs(title = "Strict arterial-only re-analysis: original EC states") + white_theme
save_plot_white(p_state, file.path(outdir, "02_plots", "Fig02_strict_arterial_only_UMAP_original_states.png"), 8, 6.5)

p_cluster <- DimPlot(obj_strict, reduction = "umap_arterial_only",
                     group.by = "strict_recluster", label = TRUE, repel = TRUE) +
  labs(title = "Strict arterial-only re-analysis: de novo clusters") + white_theme
save_plot_white(p_cluster, file.path(outdir, "02_plots", "Fig03_strict_arterial_only_UMAP_reclusters.png"), 8, 6.5)

# New cluster vs original state composition.
xt <- table(obj_strict$strict_recluster, obj_strict$arterial_state_original)
xt_count <- as.data.frame.matrix(xt)
xt_count$strict_recluster <- rownames(xt_count)
xt_count <- xt_count[, c("strict_recluster", setdiff(colnames(xt_count), "strict_recluster")), drop = FALSE]
write_tsv(xt_count, file.path(outdir, "03_tables", "06_recluster_by_original_state_counts.tsv"))

xt_prop <- prop.table(xt, margin = 1)
xt_prop_df <- as.data.frame(xt_prop)
colnames(xt_prop_df) <- c("strict_recluster", "state", "proportion")
write_tsv(xt_prop_df, file.path(outdir, "03_tables", "07_recluster_by_original_state_proportions.tsv"))

p_mix <- ggplot(xt_prop_df, aes(x = state, y = strict_recluster, fill = proportion)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", proportion)), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(x = "Original state", y = "De novo strict cluster",
       title = "Composition of de novo clusters after strict arterial-only filtering") +
  white_theme
save_plot_white(p_mix, file.path(outdir, "02_plots", "Fig04_recluster_original_state_composition.png"), 8, 5.5)

# =========================
# 6. Differential expression: retained C3 vs C7 and C3 vs C5
# =========================
Idents(obj_strict) <- factor(obj_strict$arterial_state_original, levels = state_factor_levels)

run_de <- function(a, b, outfile) {
  na <- sum(Idents(obj_strict) == a, na.rm = TRUE)
  nb <- sum(Idents(obj_strict) == b, na.rm = TRUE)
  if (na < 20 || nb < 20) {
    cat("Skipping DE", a, "vs", b, "because cell count is too small:", na, nb, "\n")
    return(NULL)
  }
  de <- FindMarkers(
    obj_strict,
    ident.1 = a,
    ident.2 = b,
    assay = "RNA",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.05,
    only.pos = FALSE,
    verbose = FALSE
  )
  de$gene <- rownames(de)
  de <- de[, c("gene", setdiff(colnames(de), "gene")), drop = FALSE]
  write.csv(de, outfile, row.names = FALSE)
  de
}

de_c3_c7 <- run_de("C3", "C7", file.path(outdir, "03_tables", "08_C3_vs_C7_DEG_strict_arterial_only.csv"))
de_c3_c5 <- run_de("C3", "C5", file.path(outdir, "03_tables", "09_C3_vs_C5_DEG_strict_arterial_only.csv"))

# Selected biological anchors from the project for easy inspection.
selected_genes <- c(
  "GJA4", "GJA5", "SEMA3G", "EFNB2", "SOX17", "HEY1", "CXCL12",
  "ANGPT2", "PGF", "KLF2", "KLF4", "THBD", "EPAS1",
  "SULF1", "ITLN1", "MGP", "ELN", "COL8A1",
  "SPARC", "FN1", "COL6A2", "IGFBP5",
  "ACKR1", "NR2F2", "RGCC", "CA4"
)
selected_genes2 <- match_genes_ci(selected_genes, rownames(obj_strict))
if (length(selected_genes2) > 0) {
  expr <- fetch_data_layer(obj_strict, selected_genes2)
  tmp <- data.frame(
    state = as.character(obj_strict@meta.data[rownames(expr), "arterial_state_original"]),
    expr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out_list <- list()
  ii <- 1
  for (st in state_factor_levels) {
    ids <- which(tmp$state == st)
    if (length(ids) == 0) next
    for (g in selected_genes2) {
      v <- as.numeric(tmp[ids, g])
      out_list[[ii]] <- data.frame(
        state = st, gene = g, n_cells = length(v),
        mean_logexpr = mean(v, na.rm = TRUE),
        median_logexpr = median(v, na.rm = TRUE),
        pct_expressing = mean(v > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      ii <- ii + 1
    }
  }
  selected_summary <- do.call(rbind, out_list)
  write_tsv(selected_summary, file.path(outdir, "03_tables", "10_selected_gene_expression_by_state.tsv"))
}

# =========================
# 7. Slingshot trajectory on strict arterial-only cells
# =========================
emb <- Embeddings(obj_strict, reduction = "umap_arterial_only")
cl <- factor(obj_strict$arterial_state_original, levels = state_factor_levels)
valid <- !is.na(cl)

# Initialize outputs so the script still finishes even if a terminal state is too sparse
# or Slingshot cannot fit a curve after the stringent filter.
obj_strict$strict_pseudotime <- rep(NA_real_, ncol(obj_strict))
trajectory_status <- "not_run"
sds <- NULL
lin_idx <- NA_integer_

state_n_strict <- table(factor(obj_strict$arterial_state_original, levels = state_factor_levels))
can_run_slingshot <- all(state_n_strict >= 10)

if (can_run_slingshot) {
  sling_result <- tryCatch({
    sds0 <- slingshot::slingshot(
      data = emb[valid, , drop = FALSE],
      clusterLabels = cl[valid],
      start.clus = "C7",
      end.clus = "C2"
    )

    ptmat <- slingshot::slingPseudotime(sds0)
    lin_list <- slingshot::slingLineages(sds0)
    writeLines(capture.output(print(lin_list)),
               file.path(outdir, "03_tables", "11_slingshot_lineages.txt"))

    if (is.null(dim(ptmat))) ptmat <- matrix(ptmat, ncol = 1)
    coverage <- apply(ptmat, 2, function(v) sum(is.finite(v)))
    li <- which.max(coverage)
    pt <- ptmat[, li]

    pt_full <- rep(NA_real_, ncol(obj_strict))
    names(pt_full) <- colnames(obj_strict)
    pt_full[rownames(emb)[valid]] <- pt
    obj_strict$strict_pseudotime <- pt_full[colnames(obj_strict)]

    list(ok = TRUE, sds = sds0, lin_idx = li)
  }, error = function(e) {
    writeLines(c("Slingshot failed after strict filtering:", conditionMessage(e)),
               file.path(outdir, "03_tables", "11_slingshot_lineages.txt"))
    list(ok = FALSE, error = conditionMessage(e))
  })

  if (isTRUE(sling_result$ok)) {
    trajectory_status <- "success"
    sds <- sling_result$sds
    lin_idx <- sling_result$lin_idx
  } else {
    trajectory_status <- paste0("failed: ", sling_result$error)
  }
} else {
  trajectory_status <- paste0(
    "skipped: at least one expected state has <10 retained cells; counts = ",
    paste(names(state_n_strict), as.integer(state_n_strict), sep = ":", collapse = ", ")
  )
  writeLines(trajectory_status,
             file.path(outdir, "03_tables", "11_slingshot_lineages.txt"))
}

pt_df <- data.frame(
  cell = colnames(obj_strict),
  state = obj_strict$arterial_state_original,
  pseudotime = obj_strict$strict_pseudotime,
  dataset = obj_strict$detected_dataset,
  sample = obj_strict$detected_sample,
  stringsAsFactors = FALSE
)
write_tsv(pt_df, file.path(outdir, "03_tables", "12_cell_pseudotime_strict.tsv"))

pt_state <- do.call(rbind, lapply(state_factor_levels, function(st) {
  v <- pt_df$pseudotime[pt_df$state == st & is.finite(pt_df$pseudotime)]
  data.frame(
    state = st,
    n_pseudotime = length(v),
    median_pseudotime = if (length(v) > 0) median(v) else NA_real_,
    q25 = if (length(v) > 0) as.numeric(quantile(v, 0.25)) else NA_real_,
    q75 = if (length(v) > 0) as.numeric(quantile(v, 0.75)) else NA_real_
  )
}))
write_tsv(pt_state, file.path(outdir, "03_tables", "13_pseudotime_by_state.tsv"))

if (any(is.finite(pt_df$pseudotime))) {
  # Pseudotime distribution plot.
  pt_plot_df <- pt_df[is.finite(pt_df$pseudotime) & pt_df$state %in% state_factor_levels, ]
  pt_plot_df$state <- factor(pt_plot_df$state, levels = state_factor_levels)
  p_pt_box <- ggplot(pt_plot_df, aes(x = state, y = pseudotime)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    labs(x = "EC state", y = "Slingshot pseudotime",
         title = "Pseudotime after strict arterial-only filtering") +
    white_theme
  save_plot_white(p_pt_box,
                  file.path(outdir, "02_plots", "Fig05_strict_pseudotime_by_state.png"), 8, 6)

  # UMAP with pseudotime and Slingshot curve.
  um <- as.data.frame(emb)
  colnames(um)[1:2] <- c("UMAP_1", "UMAP_2")
  um$cell <- rownames(um)
  pt_map <- obj_strict$strict_pseudotime
  names(pt_map) <- colnames(obj_strict)
  um$pseudotime <- pt_map[um$cell]

  p_pt <- ggplot(um, aes(x = UMAP_1, y = UMAP_2, colour = pseudotime)) +
    geom_point(size = 0.35, alpha = 0.8) +
    scale_colour_viridis_c(na.value = "grey85") +
    labs(title = "Strict arterial-only Slingshot trajectory", colour = "Pseudotime") +
    white_theme

  if (!is.null(sds) && is.finite(lin_idx)) {
    curves <- slingshot::slingCurves(sds)
    if (length(curves) >= lin_idx && !is.null(curves[[lin_idx]]$s)) {
      cv <- as.data.frame(curves[[lin_idx]]$s)
      if (ncol(cv) >= 2) {
        colnames(cv)[1:2] <- c("UMAP_1", "UMAP_2")
        p_pt <- p_pt + geom_path(
          data = cv, aes(x = UMAP_1, y = UMAP_2),
          inherit.aes = FALSE, linewidth = 1.0, colour = "black"
        )
      }
    }
  }
  save_plot_white(p_pt,
                  file.path(outdir, "02_plots", "Fig06_strict_arterial_only_slingshot.png"),
                  8, 6.5)
}

# Save final object with pseudotime (or NA pseudotime if Slingshot could not be fit).
saveRDS(obj_strict,
        file.path(outdir, "01_objects", "AS_EC_strict_arterial_only_reclustered_slingshot.rds"),
        compress = FALSE)

# =========================
# 8. Automated summary metrics
# =========================
get_count_metric <- function(filter_name, state, field) {
  z <- state_counts[state_counts$filter == filter_name & state_counts$state == state, ]
  if (nrow(z) == 0) return(NA)
  z[[field]][1]
}

c3_before <- get_count_metric("primary_arterial_only", "C3", "n_before")
c3_after <- get_count_metric("primary_arterial_only", "C3", "n_after")
c3_ret <- get_count_metric("primary_arterial_only", "C3", "retention_fraction")
c3_ultra_after <- get_count_metric("ultra_arterial_dominant", "C3", "n_after")
c3_ultra_ret <- get_count_metric("ultra_arterial_dominant", "C3", "retention_fraction")

strict_c3_md <- obj_strict@meta.data[obj_strict$arterial_state_original == "C3", , drop = FALSE]
n_c3_datasets <- length(unique(strict_c3_md$detected_dataset))
n_c3_samples <- length(unique(strict_c3_md$detected_sample))

pt_named <- setNames(pt_state$median_pseudotime, pt_state$state)
full_order_ok <- FALSE
c3_between_ok <- FALSE
if (all(state_factor_levels %in% names(pt_named)) && all(is.finite(pt_named[state_factor_levels]))) {
  full_order_ok <- all(diff(pt_named[state_factor_levels]) > 0)
  c3_between_ok <- (pt_named["C7"] < pt_named["C3"]) && (pt_named["C3"] < pt_named["C5"])
}
if (any(is.finite(pt_state$median_pseudotime))) {
  observed_order <- paste(pt_state$state[order(pt_state$median_pseudotime, na.last = NA)], collapse = " -> ")
} else {
  observed_order <- "not available"
}

summary_lines <- c(
  "STRICT ARTERIAL-ONLY C3 SENSITIVITY ANALYSIS",
  paste0("Input: ", input_rds),
  paste0("Detected state column: ", state_col),
  paste0("Detected dataset column: ", ifelse(is.null(dataset_col), "not found; sample used", dataset_col)),
  paste0("Detected sample column: ", ifelse(is.null(sample_col), "not found", sample_col)),
  "",
  "Primary filter definition:",
  paste0("  Original state in ", paste(axis_states, collapse = "/")),
  paste0("  venous_program_z <= C7 95th percentile = ", signif(ven_cut, 5)),
  paste0("  capillary_program_z <= C7 95th percentile = ", signif(cap_cut, 5)),
  "Secondary ultra-strict definition:",
  "  Primary filter + arterial_program_z > venous_program_z and > capillary_program_z",
  "",
  paste0("C3 before filtering: ", c3_before),
  paste0("C3 after primary arterial-only filtering: ", c3_after,
         " (retention ", sprintf("%.1f%%", 100 * c3_ret), ")"),
  paste0("C3 after ultra arterial-dominant filtering: ", c3_ultra_after,
         " (retention ", sprintf("%.1f%%", 100 * c3_ultra_ret), ")"),
  paste0("Retained C3 spans ", n_c3_datasets, " datasets and ", n_c3_samples, " samples."),
  "",
  paste0("Slingshot status: ", trajectory_status),
  paste0("Observed median pseudotime order: ", observed_order),
  paste0("C7 < C3 < C5 by median pseudotime: ", c3_between_ok),
  paste0("Exact C7 < C3 < C5 < C1 < C2 median order: ", full_order_ok),
  "",
  "Interpretation rule:",
  "  Strong support against venous/capillary contamination requires C3 to remain after the primary filter,",
  "  to occur across multiple independent samples/datasets, to retain its molecular contrast versus C7/C5,",
  "  and to remain positioned between C7 and C5 in the strict arterial-only trajectory.",
  "  The ultra-strict arterial-dominant result is a secondary stress test and should not be used alone,",
  "  because dysfunctional arterial cells may physiologically lose arterial-identity gene expression.",
  "",
  paste0("End time: ", format(Sys.time()))
)
writeLines(summary_lines, file.path(outdir, "03_tables", "14_arterial_only_C3_sensitivity_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# Session info for reproducibility.
writeLines(capture.output(sessionInfo()), file.path(outdir, "04_logs", "sessionInfo.txt"))

cat("\nAnalysis completed successfully.\n")
cat("Main result directory:", outdir, "\n")
