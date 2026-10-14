suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260802)

base <- "/root/GSE247238_external_validation2/09_EC_identity_contamination_sensitivity"
pb_root <- file.path(base, "02_sensitivity_pseudobulk")
statdir <- file.path(base, "03_statistics")
plotdir <- file.path(base, "plots")
dir.create(statdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

sets <- c("S1_HM_ALL", "S2_HIGH", "S3_STRICT_EC")
min_cells_by_set <- c(S1_HM_ALL = 10, S2_HIGH = 10, S3_STRICT_EC = 5)
min_paired_patients <- 3
primary_patients <- paste0("CAR", 2:10)

states <- c("C7_quiescent", "C3_early_activation", "C5_transition", "C1C2_late_remodeling")
contrasts <- list(
  C3_vs_C7 = c("C3_early_activation", "C7_quiescent"),
  Late_vs_C3 = c("C1C2_late_remodeling", "C3_early_activation"),
  Late_vs_C7 = c("C1C2_late_remodeling", "C7_quiescent")
)

candidates <- c(
  "TGFBR3", "SMAD1", "SOX4", "HIF1A", "ANGPT2",
  "NAMPT", "THBS1", "CTHRC1", "LOX"
)
target_genes <- c(candidates, "SULF1")
core_genes <- c("ANGPT2", "LOX", "CTHRC1", "SULF1")

read_gene_vector <- function(path) {
  scan(path, what = character(), sep = "\n", quiet = TRUE)
}

read_pb <- function(set_name) {
  d <- file.path(pb_root, set_name)
  x <- Matrix::readMM(file.path(d, "01_patient_state_counts.mtx"))
  x <- as(x, "dgCMatrix")
  genes <- read_gene_vector(file.path(d, "02_genes.txt"))
  meta <- read.csv(file.path(d, "03_metadata.csv"), check.names = FALSE)
  if (nrow(x) != length(genes)) stop(set_name, ": matrix rows != genes")
  if (ncol(x) != nrow(meta)) stop(set_name, ": matrix columns != metadata rows")
  rownames(x) <- toupper(genes)
  colnames(x) <- meta$pb_sample_id
  list(counts = x, meta = meta)
}

run_paired_contrast <- function(counts, meta, set_name, target_state, reference_state, contrast_name) {
  min_cells <- unname(min_cells_by_set[set_name])
  eligible <- meta[
    meta$mapped_state %in% c(reference_state, target_state) &
      meta$n_cells >= min_cells &
      meta$patient %in% primary_patients,
    , drop = FALSE
  ]
  ref_patients <- eligible$patient[eligible$mapped_state == reference_state]
  tar_patients <- eligible$patient[eligible$mapped_state == target_state]
  paired <- sort(intersect(ref_patients, tar_patients))

  if (length(paired) < min_paired_patients) {
    return(list(
      result = NULL,
      n = length(paired),
      patients = paired,
      status = "paired patients insufficient",
      min_cells = min_cells
    ))
  }

  selected <- do.call(rbind, lapply(paired, function(p) {
    rbind(
      eligible[eligible$patient == p & eligible$mapped_state == reference_state, , drop = FALSE][1, ],
      eligible[eligible$patient == p & eligible$mapped_state == target_state, , drop = FALSE][1, ]
    )
  }))
  selected$patient_factor <- factor(selected$patient)
  selected$state_factor <- factor(selected$mapped_state, levels = c(reference_state, target_state))
  x <- counts[, selected$pb_sample_id, drop = FALSE]
  design <- model.matrix(~ patient_factor + state_factor, data = selected)
  coef_idx <- grep("^state_factor", colnames(design))
  if (length(coef_idx) != 1) stop("Cannot identify state coefficient: ", set_name, " ", contrast_name)

  y <- edgeR::DGEList(counts = x)
  keep <- edgeR::filterByExpr(y, design = design)
  if (sum(keep) < 50) {
    return(list(
      result = NULL,
      n = length(paired),
      patients = paired,
      status = "too few genes after filterByExpr",
      min_cells = min_cells
    ))
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$sensitivity_set <- set_name
  tab$contrast <- contrast_name
  tab$target_state <- target_state
  tab$reference_state <- reference_state
  tab$n_paired_patients <- length(paired)
  tab$paired_patients <- paste(paired, collapse = ";")
  tab$minimum_cells_per_patient_state <- min_cells
  tab <- tab[, c(
    "gene", "sensitivity_set", "contrast", "target_state", "reference_state",
    "n_paired_patients", "paired_patients", "minimum_cells_per_patient_state",
    setdiff(colnames(tab), c(
      "gene", "sensitivity_set", "contrast", "target_state", "reference_state",
      "n_paired_patients", "paired_patients", "minimum_cells_per_patient_state"
    ))
  )]
  list(result = tab, n = length(paired), patients = paired, status = "ok", min_cells = min_cells)
}

all_results <- list()
qc_rows <- list()
for (set_name in sets) {
  pb <- read_pb(set_name)
  for (contrast_name in names(contrasts)) {
    pair <- contrasts[[contrast_name]]
    res <- tryCatch(
      run_paired_contrast(pb$counts, pb$meta, set_name, pair[1], pair[2], contrast_name),
      error = function(e) list(
        result = NULL, n = 0, patients = character(),
        status = paste0("ERROR: ", e$message),
        min_cells = unname(min_cells_by_set[set_name])
      )
    )
    qc_rows[[paste(set_name, contrast_name, sep = "|||")]] <- data.frame(
      sensitivity_set = set_name,
      contrast = contrast_name,
      n_paired_patients = res$n,
      paired_patients = paste(res$patients, collapse = ";"),
      minimum_cells_per_patient_state = res$min_cells,
      status = res$status
    )
    if (!is.null(res$result)) {
      all_results[[paste(set_name, contrast_name, sep = "|||")]] <- res$result
      write.csv(
        res$result,
        file.path(statdir, paste0("edgeR_", set_name, "_", contrast_name, "_full.csv")),
        row.names = FALSE
      )
    }
  }
}

qc <- do.call(rbind, qc_rows)
rownames(qc) <- NULL
write.csv(qc, file.path(statdir, "01_state_contrast_QC_by_sensitivity_set.csv"), row.names = FALSE)

if (length(all_results) > 0) {
  all_tab <- do.call(rbind, all_results)
  target_tab <- all_tab[all_tab$gene %in% target_genes, , drop = FALSE]
} else {
  all_tab <- data.frame()
  target_tab <- data.frame()
}
write.csv(target_tab, file.path(statdir, "02_candidate_SULF1_sensitivity_edgeR_long.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# 患者级logCPM配对差值，直接展示每位患者方向。
# --------------------------------------------------------------------------
patient_diff_rows <- list()
direction_rows <- list()
for (set_name in sets) {
  expr <- read.csv(file.path(pb_root, set_name, "04_target_gene_logCPM.csv"), check.names = FALSE)
  min_cells <- unname(min_cells_by_set[set_name])
  expr <- expr[
    expr$patient %in% primary_patients &
      expr$n_cells >= min_cells &
      expr$mapped_state %in% states,
    , drop = FALSE
  ]
  for (contrast_name in names(contrasts)) {
    target_state <- contrasts[[contrast_name]][1]
    reference_state <- contrasts[[contrast_name]][2]
    for (gene in target_genes) {
      z <- expr[expr$gene == gene & expr$mapped_state %in% c(reference_state, target_state), ]
      wide <- reshape(
        z[, c("patient", "mapped_state", "log2CPM_plus1")],
        idvar = "patient", timevar = "mapped_state", direction = "wide"
      )
      ref_col <- paste0("log2CPM_plus1.", reference_state)
      tar_col <- paste0("log2CPM_plus1.", target_state)
      if (!(ref_col %in% colnames(wide)) || !(tar_col %in% colnames(wide))) next
      wide <- wide[is.finite(wide[[ref_col]]) & is.finite(wide[[tar_col]]), , drop = FALSE]
      if (nrow(wide) == 0) next
      wide$difference_target_minus_reference <- wide[[tar_col]] - wide[[ref_col]]
      for (i in seq_len(nrow(wide))) {
        patient_diff_rows[[length(patient_diff_rows) + 1]] <- data.frame(
          sensitivity_set = set_name,
          contrast = contrast_name,
          gene = gene,
          patient = wide$patient[i],
          reference_log2CPM_plus1 = wide[[ref_col]][i],
          target_log2CPM_plus1 = wide[[tar_col]][i],
          difference_target_minus_reference = wide$difference_target_minus_reference[i]
        )
      }
      d <- wide$difference_target_minus_reference
      direction_rows[[length(direction_rows) + 1]] <- data.frame(
        sensitivity_set = set_name,
        contrast = contrast_name,
        gene = gene,
        n_evaluable_patients = length(d),
        n_positive = sum(d > 0),
        n_negative = sum(d < 0),
        n_zero = sum(d == 0),
        fraction_positive = mean(d > 0),
        median_difference = median(d)
      )
    }
  }
}

patient_diffs <- if (length(patient_diff_rows)) do.call(rbind, patient_diff_rows) else data.frame()
direction_summary <- if (length(direction_rows)) do.call(rbind, direction_rows) else data.frame()
write.csv(patient_diffs, file.path(statdir, "03_patient_target_gene_state_differences.csv"), row.names = FALSE)
write.csv(direction_summary, file.path(statdir, "04_patient_direction_summary_by_sensitivity_set.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# 核心候选稳健性判定：强调方向保留，不把FDR作为唯一标准。
# --------------------------------------------------------------------------
get_edge_value <- function(gene, set_name, contrast_name, field) {
  if (nrow(target_tab) == 0) return(NA_real_)
  z <- target_tab[
    target_tab$gene == gene &
      target_tab$sensitivity_set == set_name &
      target_tab$contrast == contrast_name,
    , drop = FALSE
  ]
  if (nrow(z) == 0 || !(field %in% colnames(z))) return(NA_real_)
  z[[field]][1]
}

robust_rows <- list()
for (gene in target_genes) {
  c3_fc <- sapply(sets, function(s) get_edge_value(gene, s, "C3_vs_C7", "logFC"))
  late_c3_fc <- sapply(sets, function(s) get_edge_value(gene, s, "Late_vs_C3", "logFC"))
  c3_available <- is.finite(c3_fc)
  late_available <- is.finite(late_c3_fc)
  n_c3_pos <- sum(c3_fc[c3_available] > 0)
  n_late_c3_pos <- sum(late_c3_fc[late_available] > 0)
  n_late_c3_neg <- sum(late_c3_fc[late_available] < 0)

  strict_c3 <- unname(c3_fc["S3_STRICT_EC"])
  strict_late <- unname(late_c3_fc["S3_STRICT_EC"])

  classification <- if (gene == "SULF1") {
    if (sum(late_available) >= 2 && n_late_c3_pos == sum(late_available)) {
      "late_anchor_direction_preserved"
    } else if (sum(late_available) >= 1 && n_late_c3_pos >= 1) {
      "late_anchor_partial_support"
    } else {
      "late_anchor_not_supported_or_insufficient"
    }
  } else if (sum(c3_available) >= 2 && n_c3_pos == sum(c3_available)) {
    if (sum(late_available) >= 2 && n_late_c3_neg == sum(late_available)) {
      "early_C3_peak_direction_preserved"
    } else {
      "early_C3_activation_direction_preserved"
    }
  } else if (sum(c3_available) >= 1 && n_c3_pos >= 1) {
    "partial_or_set_dependent_support"
  } else {
    "not_supported_or_insufficient"
  }

  robust_rows[[gene]] <- data.frame(
    gene = gene,
    n_sets_C3_vs_C7_available = sum(c3_available),
    n_sets_C3_vs_C7_positive = n_c3_pos,
    n_sets_Late_vs_C3_available = sum(late_available),
    n_sets_Late_vs_C3_positive = n_late_c3_pos,
    n_sets_Late_vs_C3_negative = n_late_c3_neg,
    strict_EC_C3_vs_C7_logFC = strict_c3,
    strict_EC_Late_vs_C3_logFC = strict_late,
    direction_classification = classification,
    interpretation = "patient-level mapped-state sensitivity; not longitudinal timing or causal proof"
  )
}
robust <- do.call(rbind, robust_rows)
rownames(robust) <- NULL
write.csv(robust, file.path(statdir, "05_candidate_robustness_summary.csv"), row.names = FALSE)
write.csv(
  robust[robust$gene %in% core_genes, , drop = FALSE],
  file.path(statdir, "06_core_ANGPT2_LOX_CTHRC1_SULF1_summary.csv"),
  row.names = FALSE
)

# --------------------------------------------------------------------------
# 纯白背景图。
# --------------------------------------------------------------------------
white_theme <- theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

if (nrow(target_tab) > 0) {
  heat <- target_tab[target_tab$contrast %in% c("C3_vs_C7", "Late_vs_C3"), ]
  heat$gene <- factor(heat$gene, levels = target_genes)
  heat$set_contrast <- paste(heat$sensitivity_set, heat$contrast, sep = " | ")
  p <- ggplot(heat, aes(x = set_contrast, y = gene, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", logFC)), size = 2.8) +
    scale_fill_gradient2(midpoint = 0) +
    labs(
      title = "Candidate state logFC across EC-identity sensitivity filters",
      subtitle = "Patient-paired edgeR; positive C3_vs_C7 supports early activation; positive Late_vs_C3 supports late enrichment",
      x = NULL, y = NULL, fill = "logFC"
    ) + white_theme
  ggsave(file.path(plotdir, "06_candidate_edgeR_logFC_sensitivity_heatmap.png"), p,
         width = 13, height = 7.5, dpi = 300, bg = "white")
}

if (nrow(patient_diffs) > 0) {
  z <- patient_diffs[
    patient_diffs$gene %in% core_genes &
      patient_diffs$contrast %in% c("C3_vs_C7", "Late_vs_C3"),
    , drop = FALSE
  ]
  z$gene <- factor(z$gene, levels = core_genes)
  p <- ggplot(z, aes(x = sensitivity_set, y = difference_target_minus_reference)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.8, alpha = 0.8) +
    stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.4) +
    facet_grid(gene ~ contrast, scales = "free_y") +
    labs(
      title = "Patient-level core-gene state differences across filters",
      x = "Sensitivity set", y = "Target state - reference state log2(CPM+1)"
    ) + white_theme
  ggsave(file.path(plotdir, "07_core_gene_patient_differences_across_filters.png"), p,
         width = 11, height = 10, dpi = 300, bg = "white")
}

qc_text <- paste(capture.output(print(qc, row.names = FALSE)), collapse = "\n")
core_text <- paste(capture.output(print(robust[robust$gene %in% core_genes, ], row.names = FALSE)), collapse = "\n")
summary_lines <- c(
  "GSE247238 Step 5 sensitivity edgeR completed.",
  "",
  "Statistical unit: patient.",
  "S1_HM_ALL and S2_HIGH require >=10 cells per patient-state; S3_STRICT_EC requires >=5 cells per patient-state.",
  "C5 was not tested because the external dataset contained too few C5-like cells.",
  "",
  "Contrast availability:",
  qc_text,
  "",
  "Core-gene direction summary:",
  core_text,
  "",
  "Interpretation limits:",
  "1. Direction preserved after strict filtering argues against a signal driven solely by obvious non-EC co-expression.",
  "2. Loss of significance after filtering may reflect lower cell number as well as contamination; inspect effect direction and patient consistency.",
  "3. This analysis does not prove lineage, longitudinal timing, or causality.",
  "4. Non-EC risk modules are QC indicators, not direct doublet calls."
)
writeLines(summary_lines, file.path(statdir, "README_statistics_and_interpretation.txt"))
writeLines(capture.output(sessionInfo()), file.path(statdir, "sessionInfo_R.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
