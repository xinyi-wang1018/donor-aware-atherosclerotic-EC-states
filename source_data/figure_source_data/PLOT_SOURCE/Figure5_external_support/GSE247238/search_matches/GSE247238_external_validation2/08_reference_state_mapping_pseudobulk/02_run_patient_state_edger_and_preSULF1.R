suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")
set.seed(20260801)

base <- "/root/GSE247238_external_validation2/08_reference_state_mapping_pseudobulk"
pbdir <- file.path(base, "03_patient_pseudobulk")
plotdir <- file.path(base, "plots")
dir.create(pbdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS_PRIMARY <- 10
MIN_PAIRED_PATIENTS <- 3

trajectory_states <- c(
  "C7_quiescent", "C3_early_activation", "C5_transition", "C1C2_late_remodeling"
)
state_order <- setNames(seq_along(trajectory_states), trajectory_states)

candidates <- c(
  "TGFBR3", "SMAD1", "SOX4", "HIF1A", "ANGPT2",
  "NAMPT", "THBS1", "CTHRC1", "LOX"
)
target_genes <- c(candidates, "SULF1")

read_gene_vector <- function(path) {
  scan(path, what = character(), sep = "\n", quiet = TRUE)
}

read_pb <- function(mtx_path, gene_path, meta_path) {
  x <- Matrix::readMM(mtx_path)
  x <- as(x, "dgCMatrix")
  genes <- read_gene_vector(gene_path)
  meta <- read.csv(meta_path, check.names = FALSE)
  if (nrow(x) != length(genes)) stop("pseudobulk矩阵行数与基因数不一致。")
  if (ncol(x) != nrow(meta)) stop("pseudobulk矩阵列数与metadata不一致。")
  rownames(x) <- genes
  colnames(x) <- meta$pb_sample_id
  list(counts = x, meta = meta)
}

ps <- read_pb(
  file.path(pbdir, "01_patient_state_counts.mtx"),
  file.path(pbdir, "02_patient_state_genes.txt"),
  file.path(pbdir, "03_patient_state_metadata.csv")
)
pcs <- read_pb(
  file.path(pbdir, "04_patient_condition_state_counts.mtx"),
  file.path(pbdir, "05_patient_condition_state_genes.txt"),
  file.path(pbdir, "06_patient_condition_state_metadata.csv")
)

# 大小写统一并去除重复基因；重复时求和。
collapse_duplicate_genes <- function(counts) {
  upper <- toupper(rownames(counts))
  if (!anyDuplicated(upper)) {
    rownames(counts) <- upper
    return(counts)
  }
  f <- factor(upper, levels = unique(upper))
  agg <- Matrix::sparse.model.matrix(~ 0 + f)
  colnames(agg) <- levels(f)
  out <- Matrix::t(agg) %*% counts
  rownames(out) <- levels(f)
  out
}
ps$counts <- collapse_duplicate_genes(ps$counts)
pcs$counts <- collapse_duplicate_genes(pcs$counts)

run_paired_state_contrast <- function(counts, meta, target_state, reference_state, contrast_name) {
  eligible <- meta[
    meta$mapped_state %in% c(reference_state, target_state) &
      meta$n_cells >= MIN_CELLS_PRIMARY &
      meta$patient %in% paste0("CAR", 2:10),
    , drop = FALSE
  ]
  ref_patients <- eligible$patient[eligible$mapped_state == reference_state]
  tar_patients <- eligible$patient[eligible$mapped_state == target_state]
  paired <- sort(intersect(ref_patients, tar_patients))

  if (length(paired) < MIN_PAIRED_PATIENTS) {
    return(list(result = NULL, meta = NULL, n = length(paired), message = "paired patients insufficient"))
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
  if (length(coef_idx) != 1) stop("无法确定状态对比系数：", contrast_name)

  y <- edgeR::DGEList(counts = x)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$contrast <- contrast_name
  tab$target_state <- target_state
  tab$reference_state <- reference_state
  tab$n_paired_patients <- length(paired)
  tab$paired_patients <- paste(paired, collapse = ";")
  tab <- tab[, c("gene", "contrast", "target_state", "reference_state", "n_paired_patients", "paired_patients", setdiff(colnames(tab), c("gene", "contrast", "target_state", "reference_state", "n_paired_patients", "paired_patients")))]
  list(result = tab, meta = selected, n = length(paired), message = "ok")
}

state_contrasts <- list(
  C3_vs_C7 = c("C3_early_activation", "C7_quiescent"),
  C5_vs_C3 = c("C5_transition", "C3_early_activation"),
  C5_vs_C7 = c("C5_transition", "C7_quiescent"),
  Late_vs_C5 = c("C1C2_late_remodeling", "C5_transition"),
  Late_vs_C3 = c("C1C2_late_remodeling", "C3_early_activation"),
  Late_vs_C7 = c("C1C2_late_remodeling", "C7_quiescent")
)

state_result_list <- list()
state_qc <- list()
for (nm in names(state_contrasts)) {
  pair <- state_contrasts[[nm]]
  res <- tryCatch(
    run_paired_state_contrast(ps$counts, ps$meta, pair[1], pair[2], nm),
    error = function(e) list(result = NULL, meta = NULL, n = 0, message = paste("ERROR:", e$message))
  )
  state_qc[[nm]] <- data.frame(contrast = nm, n_paired_patients = res$n, status = res$message)
  if (!is.null(res$result)) {
    state_result_list[[nm]] <- res$result
    write.csv(res$result, file.path(pbdir, paste0("10_edgeR_", nm, "_full.csv")), row.names = FALSE)
  }
}
state_qc_df <- do.call(rbind, state_qc)
write.csv(state_qc_df, file.path(pbdir, "10_state_contrast_QC.csv"), row.names = FALSE)

if (length(state_result_list) > 0) {
  all_state_results <- do.call(rbind, state_result_list)
  target_state_results <- all_state_results[all_state_results$gene %in% target_genes, , drop = FALSE]
  write.csv(target_state_results, file.path(pbdir, "11_candidate_SULF1_state_contrast_edgeR.csv"), row.names = FALSE)
} else {
  all_state_results <- data.frame()
  target_state_results <- data.frame()
  write.csv(target_state_results, file.path(pbdir, "11_candidate_SULF1_state_contrast_edgeR.csv"), row.names = FALSE)
}

# 同一映射状态内，患者配对Plaque-Control edgeR。
run_condition_within_state <- function(counts, meta, state_name) {
  eligible <- meta[
    meta$mapped_state == state_name &
      meta$n_cells >= MIN_CELLS_PRIMARY &
      meta$patient %in% paste0("CAR", 2:10) &
      meta$condition %in% c("Control", "Plaque"),
    , drop = FALSE
  ]
  ctrl <- eligible$patient[eligible$condition == "Control"]
  plaq <- eligible$patient[eligible$condition == "Plaque"]
  paired <- sort(intersect(ctrl, plaq))
  if (length(paired) < MIN_PAIRED_PATIENTS) {
    return(list(result = NULL, n = length(paired), message = "paired patients insufficient"))
  }
  selected <- do.call(rbind, lapply(paired, function(p) {
    rbind(
      eligible[eligible$patient == p & eligible$condition == "Control", , drop = FALSE][1, ],
      eligible[eligible$patient == p & eligible$condition == "Plaque", , drop = FALSE][1, ]
    )
  }))
  selected$patient_factor <- factor(selected$patient)
  selected$condition_factor <- factor(selected$condition, levels = c("Control", "Plaque"))
  x <- counts[, selected$pb_sample_id, drop = FALSE]
  design <- model.matrix(~ patient_factor + condition_factor, data = selected)
  coef_idx <- grep("^condition_factor", colnames(design))
  y <- edgeR::DGEList(counts = x)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$mapped_state <- state_name
  tab$contrast <- "Plaque_vs_Control_within_state"
  tab$n_paired_patients <- length(paired)
  tab$paired_patients <- paste(paired, collapse = ";")
  list(result = tab, n = length(paired), message = "ok")
}

condition_results <- list()
condition_qc <- list()
for (s in trajectory_states) {
  res <- tryCatch(
    run_condition_within_state(pcs$counts, pcs$meta, s),
    error = function(e) list(result = NULL, n = 0, message = paste("ERROR:", e$message))
  )
  condition_qc[[s]] <- data.frame(mapped_state = s, n_paired_patients = res$n, status = res$message)
  if (!is.null(res$result)) condition_results[[s]] <- res$result
}
write.csv(do.call(rbind, condition_qc), file.path(pbdir, "12_within_state_condition_QC.csv"), row.names = FALSE)
if (length(condition_results) > 0) {
  condition_all <- do.call(rbind, condition_results)
  condition_target <- condition_all[condition_all$gene %in% target_genes, , drop = FALSE]
  write.csv(condition_target, file.path(pbdir, "13_candidate_SULF1_Plaque_vs_Control_within_state_edgeR.csv"), row.names = FALSE)
} else {
  condition_target <- data.frame()
  write.csv(condition_target, file.path(pbdir, "13_candidate_SULF1_Plaque_vs_Control_within_state_edgeR.csv"), row.names = FALSE)
}

# --------------------------------------------------------------------------
# 描述性“早于SULF1”状态顺序判定。
# 这是状态顺序支持，不是时间因果或真实纵向追踪。
# --------------------------------------------------------------------------
expr <- read.csv(file.path(pbdir, "07_patient_state_candidate_SULF1_logCPM.csv"), check.names = FALSE)
expr <- expr[
  expr$patient %in% paste0("CAR", 2:10) &
    expr$mapped_state %in% trajectory_states &
    expr$n_cells >= MIN_CELLS_PRIMARY,
  , drop = FALSE
]
expr$state_order <- state_order[expr$mapped_state]

# 固定统一规则：
# 1) CPM >= 1；
# 2) logCPM达到该患者该基因动态范围的35%；
# 3) 且至少比C7高0.25 log2(CPM+1)（若C7可用）。
get_onset <- function(df) {
  df <- df[order(df$state_order), , drop = FALSE]
  if (nrow(df) < 2 || all(!is.finite(df$log2CPM_plus1))) {
    return(data.frame(onset_state = NA_character_, onset_order = NA_real_, n_states = nrow(df)))
  }
  vals <- df$log2CPM_plus1
  vmin <- min(vals, na.rm = TRUE)
  vmax <- max(vals, na.rm = TRUE)
  threshold_dynamic <- vmin + 0.35 * (vmax - vmin)
  c7_value <- df$log2CPM_plus1[df$mapped_state == "C7_quiescent"]
  c7_threshold <- if (length(c7_value) && is.finite(c7_value[1])) c7_value[1] + 0.25 else -Inf
  threshold <- max(threshold_dynamic, c7_threshold)
  eligible <- which(df$CPM >= 1 & df$log2CPM_plus1 >= threshold)
  if (!length(eligible)) {
    return(data.frame(onset_state = NA_character_, onset_order = NA_real_, n_states = nrow(df)))
  }
  i <- eligible[1]
  data.frame(onset_state = df$mapped_state[i], onset_order = df$state_order[i], n_states = nrow(df))
}

onset_rows <- list()
for (p in unique(expr$patient)) {
  for (g in target_genes) {
    z <- expr[expr$patient == p & expr$gene == g, , drop = FALSE]
    onset <- get_onset(z)
    onset_rows[[paste(p, g, sep = "|||")]] <- data.frame(
      patient = p,
      gene = g,
      onset_state = onset$onset_state,
      onset_order = onset$onset_order,
      n_states_evaluable = onset$n_states,
      stringsAsFactors = FALSE
    )
  }
}
onset_df <- do.call(rbind, onset_rows)
rownames(onset_df) <- NULL
write.csv(onset_df, file.path(pbdir, "14_patient_gene_state_onset_descriptive.csv"), row.names = FALSE)

sulf1_onset <- onset_df[onset_df$gene == "SULF1", c("patient", "onset_state", "onset_order")]
colnames(sulf1_onset) <- c("patient", "SULF1_onset_state", "SULF1_onset_order")
compare <- merge(onset_df[onset_df$gene %in% candidates, ], sulf1_onset, by = "patient", all.x = TRUE)
compare$evaluable_vs_SULF1 <- is.finite(compare$onset_order) & is.finite(compare$SULF1_onset_order)
compare$earlier_than_SULF1 <- compare$evaluable_vs_SULF1 & compare$onset_order < compare$SULF1_onset_order
compare$same_as_SULF1 <- compare$evaluable_vs_SULF1 & compare$onset_order == compare$SULF1_onset_order
compare$later_than_SULF1 <- compare$evaluable_vs_SULF1 & compare$onset_order > compare$SULF1_onset_order
write.csv(compare, file.path(pbdir, "15_candidate_vs_SULF1_patient_state_precedence.csv"), row.names = FALSE)

precedence_summary <- do.call(rbind, lapply(candidates, function(g) {
  z <- compare[compare$gene == g & compare$evaluable_vs_SULF1, , drop = FALSE]
  n <- nrow(z)
  n_earlier <- sum(z$earlier_than_SULF1)
  frac <- if (n > 0) n_earlier / n else NA_real_
  support <- if (n >= 4 && frac >= 0.75) {
    "strong_state_order_support"
  } else if (n >= 3 && frac >= 0.60) {
    "moderate_state_order_support"
  } else if (n >= 2 && frac > 0.50) {
    "limited_directional_support"
  } else {
    "not_supported_or_insufficient"
  }
  data.frame(
    gene = g,
    n_evaluable_patients = n,
    n_earlier_than_SULF1 = n_earlier,
    fraction_earlier_than_SULF1 = frac,
    n_same_as_SULF1 = sum(z$same_as_SULF1),
    n_later_than_SULF1 = sum(z$later_than_SULF1),
    state_order_support = support,
    interpretation = "cross-sectional mapped-state precedence; not longitudinal timing or causality",
    stringsAsFactors = FALSE
  )
}))
write.csv(precedence_summary, file.path(pbdir, "16_candidate_preSULF1_state_order_summary.csv"), row.names = FALSE)

# 将状态edgeR结果整理成每个候选的证据矩阵。
if (nrow(target_state_results) > 0) {
  evidence <- target_state_results[, c("gene", "contrast", "logFC", "PValue", "FDR", "n_paired_patients")]
  write.csv(evidence, file.path(pbdir, "17_candidate_state_edgeR_evidence_matrix_long.csv"), row.names = FALSE)
}

# 绘制患者状态pseudobulk中位表达折线图。
white_theme <- theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

plot_df <- aggregate(
  log2CPM_plus1 ~ gene + mapped_state,
  data = expr,
  FUN = median,
  na.rm = TRUE
)
plot_df$mapped_state <- factor(plot_df$mapped_state, levels = trajectory_states)
plot_df$gene <- factor(plot_df$gene, levels = target_genes)

p <- ggplot(plot_df, aes(x = mapped_state, y = log2CPM_plus1, group = gene)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  labs(
    title = "GSE247238 patient-state pseudobulk profiles",
    subtitle = "Median across patients with at least 10 mapped cells per state",
    x = "Mapped state",
    y = "Median log2(CPM+1)"
  ) + white_theme

ggsave(
  file.path(plotdir, "05_candidate_SULF1_patient_state_profiles.png"),
  p, width = 12, height = 11, dpi = 300, bg = "white"
)

# 中文/英文摘要文件。
state_qc_text <- paste(capture.output(print(state_qc_df, row.names = FALSE)), collapse = "\n")
precedence_text <- paste(capture.output(print(precedence_summary, row.names = FALSE)), collapse = "\n")
summary_lines <- c(
  "GSE247238 Step 4 patient-level mapped-state pseudobulk completed.",
  "Primary statistical unit: patient.",
  paste0("Primary minimum cells per patient-state pseudobulk: ", MIN_CELLS_PRIMARY),
  "",
  "State edgeR contrast availability:",
  state_qc_text,
  "",
  "Candidate state-order precedence relative to SULF1:",
  precedence_text,
  "",
  "Interpretation limits:",
  "1. Mapped-state order supports transcriptional precedence but is not true longitudinal timing.",
  "2. A candidate earlier than SULF1 is not automatically a causal driver.",
  "3. C4 is excluded from the arterial trajectory and retained as an independent decoy state.",
  "4. Low-confidence cell mappings do not enter the primary pseudobulk analysis.",
  "5. edgeR results require both effect direction and patient consistency; FDR alone is not sufficient."
)
writeLines(summary_lines, file.path(pbdir, "README_patient_pseudobulk_and_preSULF1.txt"))
writeLines(capture.output(sessionInfo()), file.path(pbdir, "sessionInfo_pseudobulk_R.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
