options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 30 * 1024^3)
set.seed(20260731)

BASE <- "/root/GSE159677_9candidate_validation"
OUT_DIR <- file.path(BASE, "results")
OBJ_DIR <- file.path(BASE, "objects")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

path_lines <- readLines(file.path(BASE, "input_paths.txt"))
path_map <- setNames(sub("^[^=]+=", "", path_lines), sub("=.*$", "", path_lines))
EC_RDS <- unname(path_map[["EC_RDS"]])

required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "tidyr", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required R packages: ", paste(missing, collapse = ", "),
    "\n请进入之前可运行 Seurat 的 R 环境后重跑。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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
    filename = file.path(OUT_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
}

join_layers_safe <- function(obj, assay = "RNA") {
  tryCatch(
    SeuratObject::JoinLayers(obj, assay = assay),
    error = function(e1) tryCatch(
      Seurat::JoinLayers(obj, assay = assay),
      error = function(e2) obj
    )
  )
}

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

safe_num <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}

classify_direction <- function(n_ac, n_pa, n_tie, n_pairs) {
  if (n_pairs == 0) return("not_assessed")
  if (n_ac == n_pairs) return(paste0("AC_higher_", n_pairs, "_of_", n_pairs))
  if (n_pa == n_pairs) return(paste0("PA_higher_", n_pairs, "_of_", n_pairs))
  if (n_ac >= 2) return(paste0("AC_higher_", n_ac, "_of_", n_pairs))
  if (n_pa >= 2) return(paste0("PA_higher_", n_pa, "_of_", n_pairs))
  if (n_tie > 0) return("mixed_with_ties")
  "mixed_direction"
}

paired_stats <- function(df, value_col, gene_col = "gene", min_cells = NULL) {
  genes <- unique(df[[gene_col]])
  out <- lapply(genes, function(g) {
    dd <- df[df[[gene_col]] == g, , drop = FALSE]
    if (!is.null(min_cells) && "n_cells" %in% colnames(dd)) {
      dd <- dd[dd$n_cells >= min_cells, , drop = FALSE]
    }
    dd <- dd[, intersect(c("patient_id", "condition", value_col), colnames(dd)), drop = FALSE]
    if (!all(c("patient_id", "condition", value_col) %in% colnames(dd))) return(NULL)
    wide <- tidyr::pivot_wider(dd, names_from = condition, values_from = all_of(value_col))
    if (!all(c("PA", "AC") %in% colnames(wide))) return(NULL)
    wide <- wide[complete.cases(wide[, c("PA", "AC")]), , drop = FALSE]
    if (nrow(wide) == 0) return(NULL)
    wide <- wide[order(wide$patient_id), , drop = FALSE]
    diff <- safe_num(wide$AC - wide$PA)
    keep <- is.finite(diff)
    wide <- wide[keep, , drop = FALSE]
    diff <- diff[keep]
    if (length(diff) == 0) return(NULL)

    p <- if (length(diff) >= 2) {
      tryCatch(
        wilcox.test(wide$AC, wide$PA, paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    } else NA_real_

    dz <- if (length(diff) >= 2 && is.finite(sd(diff)) && sd(diff) > 0) {
      mean(diff) / sd(diff)
    } else NA_real_

    get_diff <- function(pid) {
      idx <- match(pid, wide$patient_id)
      if (is.na(idx)) NA_real_ else diff[idx]
    }

    n_ac <- sum(diff > 0, na.rm = TRUE)
    n_pa <- sum(diff < 0, na.rm = TRUE)
    n_tie <- sum(diff == 0, na.rm = TRUE)

    data.frame(
      gene = g,
      n_pairs = length(diff),
      mean_AC_minus_PA = mean(diff),
      median_AC_minus_PA = median(diff),
      paired_Cohen_dz = dz,
      wilcoxon_p = p,
      n_AC_higher = n_ac,
      n_PA_higher = n_pa,
      n_ties = n_tie,
      direction_pattern = classify_direction(n_ac, n_pa, n_tie, length(diff)),
      P1_AC_minus_PA = get_diff("P1"),
      P2_AC_minus_PA = get_diff("P2"),
      P3_AC_minus_PA = get_diff("P3")
    )
  })
  bind_rows(out)
}

if (!file.exists(EC_RDS)) stop("Input object not found: ", EC_RDS)
ec <- readRDS(EC_RDS)
if (!inherits(ec, "Seurat")) stop("Input object is not a Seurat object.")
if (!"RNA" %in% Assays(ec)) stop("Input object has no RNA assay.")

ec <- join_layers_safe(ec, "RNA")
DefaultAssay(ec) <- "RNA"

required_meta <- c("sample_id", "patient_id", "condition", "EC_layer_v3")
missing_meta <- setdiff(required_meta, colnames(ec@meta.data))
if (length(missing_meta) > 0) {
  stop("Input object lacks metadata columns: ", paste(missing_meta, collapse = ", "))
}

ec$sample_id <- as.character(ec$sample_id)
ec$patient_id <- as.character(ec$patient_id)
ec$condition <- factor(as.character(ec$condition), levels = c("PA", "AC"))
ec$EC_layer_v3 <- as.character(ec$EC_layer_v3)
if (any(is.na(ec$condition))) stop("condition contains values other than PA/AC.")

sample_meta <- ec@meta.data %>%
  transmute(
    sample_id = as.character(sample_id),
    patient_id = as.character(patient_id),
    condition = factor(as.character(condition), levels = c("PA", "AC"))
  ) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  arrange(patient_id, condition)

if (nrow(sample_meta) != 6) {
  warning("Expected 6 PA/AC samples, found ", nrow(sample_meta), ".")
}
write.csv(sample_meta, file.path(OUT_DIR, "00_sample_metadata.csv"), row.names = FALSE)

candidate_genes <- c(
  "TGFBR3", "SMAD1", "SOX4", "HIF1A", "ANGPT2",
  "NAMPT", "THBS1", "CTHRC1", "LOX"
)
late_anchor <- "SULF1"
target_genes <- c(candidate_genes, late_anchor)

# 大小写不敏感地匹配基因符号，但输出统一使用标准目标符号。
feature_names <- rownames(ec)
feature_lookup <- split(feature_names, toupper(feature_names))
matched_feature <- vapply(target_genes, function(g) {
  hit <- feature_lookup[[toupper(g)]]
  if (is.null(hit) || length(hit) == 0) NA_character_ else hit[1]
}, character(1))

gene_availability <- data.frame(
  gene = target_genes,
  role = ifelse(target_genes == late_anchor, "late_SULF1_anchor", "early_candidate"),
  matched_feature = unname(matched_feature),
  present_in_object = !is.na(matched_feature),
  stringsAsFactors = FALSE
)
write.csv(gene_availability, file.path(OUT_DIR, "01_candidate_gene_availability.csv"), row.names = FALSE)

present_targets <- gene_availability$gene[gene_availability$present_in_object]
if (length(present_targets) == 0) stop("None of the target genes are present in the object.")
actual_features <- setNames(gene_availability$matched_feature[gene_availability$present_in_object], present_targets)

counts <- get_assay_matrix(ec, "RNA", "counts")
data_mat <- get_assay_matrix(ec, "RNA", "data")

# =============================================================================
# A. Expanded EC：每个样本的平均表达、检出率和细胞数量
# =============================================================================
sample_expression_list <- lapply(sample_meta$sample_id, function(sid) {
  cells <- colnames(ec)[ec$sample_id == sid]
  md <- sample_meta[sample_meta$sample_id == sid, , drop = FALSE]
  if (length(cells) == 0) return(NULL)

  vals_mean <- vapply(present_targets, function(g) {
    f <- actual_features[[g]]
    as.numeric(Matrix::rowMeans(data_mat[f, cells, drop = FALSE]))
  }, numeric(1))

  vals_pct <- vapply(present_targets, function(g) {
    f <- actual_features[[g]]
    100 * as.numeric(Matrix::rowMeans(counts[f, cells, drop = FALSE] > 0))
  }, numeric(1))

  data.frame(
    sample_id = sid,
    patient_id = as.character(md$patient_id),
    condition = as.character(md$condition),
    n_cells = length(cells),
    gene = present_targets,
    mean_log_normalized_expression = vals_mean[present_targets],
    pct_cells_expressing = vals_pct[present_targets],
    stringsAsFactors = FALSE
  )
})
sample_expression <- bind_rows(sample_expression_list)
sample_expression$condition <- factor(sample_expression$condition, levels = c("PA", "AC"))
write.csv(sample_expression, file.path(OUT_DIR, "02_sample_level_expression_expanded_EC.csv"), row.names = FALSE)

avg_stats <- paired_stats(sample_expression, "mean_log_normalized_expression")
pct_stats <- paired_stats(sample_expression, "pct_cells_expressing")
write.csv(avg_stats, file.path(OUT_DIR, "03_paired_statistics_mean_expression.csv"), row.names = FALSE)
write.csv(pct_stats, file.path(OUT_DIR, "04_paired_statistics_detection_fraction.csv"), row.names = FALSE)

p_avg <- ggplot(
  sample_expression,
  aes(condition, mean_log_normalized_expression, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 2.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  labs(
    x = NULL,
    y = "Mean log-normalized expression",
    title = "Nine early candidates plus SULF1: paired PA versus AC"
  )
save_png(p_avg, "01_candidate_paired_mean_expression.png", 16, 8.5)

p_pct <- ggplot(
  sample_expression,
  aes(condition, pct_cells_expressing, group = patient_id, colour = patient_id)
) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 2.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  labs(
    x = NULL,
    y = "Cells expressing gene (%)",
    title = "Detection fraction: paired PA versus AC"
  )
save_png(p_pct, "02_candidate_paired_detection_fraction.png", 16, 8.5)

# =============================================================================
# B. Tier 1 / Tier 2：判断信号来自经典 EC 还是转化候选 EC
# =============================================================================
layer_levels <- c("Tier1_classical_EC", "Tier2_transition_candidate", "Tier2_low_marker_EC")
layer_expression_list <- list()
idx <- 1L
for (layer_name in intersect(layer_levels, unique(ec$EC_layer_v3))) {
  for (sid in sample_meta$sample_id) {
    cells <- colnames(ec)[ec$sample_id == sid & ec$EC_layer_v3 == layer_name]
    if (length(cells) == 0) next
    md <- sample_meta[sample_meta$sample_id == sid, , drop = FALSE]

    vals_mean <- vapply(present_targets, function(g) {
      f <- actual_features[[g]]
      as.numeric(Matrix::rowMeans(data_mat[f, cells, drop = FALSE]))
    }, numeric(1))

    vals_pct <- vapply(present_targets, function(g) {
      f <- actual_features[[g]]
      100 * as.numeric(Matrix::rowMeans(counts[f, cells, drop = FALSE] > 0))
    }, numeric(1))

    layer_expression_list[[idx]] <- data.frame(
      sample_id = sid,
      patient_id = as.character(md$patient_id),
      condition = as.character(md$condition),
      EC_layer_v3 = layer_name,
      n_cells = length(cells),
      gene = present_targets,
      mean_log_normalized_expression = vals_mean[present_targets],
      pct_cells_expressing = vals_pct[present_targets],
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}
layer_expression <- bind_rows(layer_expression_list)
if (nrow(layer_expression) > 0) {
  layer_expression$condition <- factor(layer_expression$condition, levels = c("PA", "AC"))
  write.csv(layer_expression, file.path(OUT_DIR, "05_sample_expression_by_EC_layer.csv"), row.names = FALSE)

  layer_avg_stats <- layer_expression %>%
    group_by(EC_layer_v3) %>%
    group_modify(~ paired_stats(.x, "mean_log_normalized_expression", min_cells = 30)) %>%
    ungroup()
  layer_pct_stats <- layer_expression %>%
    group_by(EC_layer_v3) %>%
    group_modify(~ paired_stats(.x, "pct_cells_expressing", min_cells = 30)) %>%
    ungroup()
  write.csv(layer_avg_stats, file.path(OUT_DIR, "06_paired_mean_expression_by_EC_layer_min30cells.csv"), row.names = FALSE)
  write.csv(layer_pct_stats, file.path(OUT_DIR, "07_paired_detection_by_EC_layer_min30cells.csv"), row.names = FALSE)

  p_layer <- ggplot(
    layer_expression %>% filter(EC_layer_v3 %in% c("Tier1_classical_EC", "Tier2_transition_candidate")),
    aes(condition, mean_log_normalized_expression, group = patient_id, colour = patient_id)
  ) +
    geom_line(linewidth = 0.55, alpha = 0.8) +
    geom_point(aes(size = pmax(n_cells, 1)), alpha = 0.9) +
    scale_size_continuous(range = c(1.5, 4.5)) +
    facet_grid(EC_layer_v3 ~ gene, scales = "free_y") +
    labs(
      x = NULL,
      y = "Mean log-normalized expression",
      size = "Cells",
      title = "Candidate expression in classical versus transition-candidate EC"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_png(p_layer, "03_candidate_expression_by_EC_layer.png", 22, 8.5)
}

# =============================================================================
# C. 映射状态中的表达：只用于判断位置是否相符，不代表真实时间顺序
# =============================================================================
if ("mapped_state" %in% colnames(ec@meta.data)) {
  ec$mapped_state <- as.character(ec$mapped_state)
  state_expression_list <- lapply(sort(unique(ec$mapped_state)), function(st) {
    cells <- colnames(ec)[ec$mapped_state == st]
    if (length(cells) == 0) return(NULL)
    vals_mean <- vapply(present_targets, function(g) {
      f <- actual_features[[g]]
      as.numeric(Matrix::rowMeans(data_mat[f, cells, drop = FALSE]))
    }, numeric(1))
    vals_pct <- vapply(present_targets, function(g) {
      f <- actual_features[[g]]
      100 * as.numeric(Matrix::rowMeans(counts[f, cells, drop = FALSE] > 0))
    }, numeric(1))
    data.frame(
      mapped_state = st,
      n_cells = length(cells),
      gene = present_targets,
      mean_log_normalized_expression = vals_mean[present_targets],
      pct_cells_expressing = vals_pct[present_targets],
      stringsAsFactors = FALSE
    )
  })
  state_expression <- bind_rows(state_expression_list)
  write.csv(state_expression, file.path(OUT_DIR, "08_candidate_expression_by_mapped_state.csv"), row.names = FALSE)

  dot_features <- unname(actual_features[present_targets])
  names(dot_features) <- present_targets
  p_state_dot <- DotPlot(ec, features = dot_features, group.by = "mapped_state") +
    RotatedAxis() +
    ggtitle("Nine early candidates and SULF1 across mapped EC states")
  save_png(p_state_dot, "04_candidate_dotplot_by_mapped_state.png", 15, 6.5)
}

# =============================================================================
# D. 患者级 pseudobulk：每位患者每个区域只作为一个生物学样本
# =============================================================================
samples <- sample_meta$sample_id
pb_counts <- sapply(samples, function(sid) {
  cells <- colnames(ec)[ec$sample_id == sid]
  Matrix::rowSums(counts[, cells, drop = FALSE])
})
if (!is.matrix(pb_counts)) pb_counts <- as.matrix(pb_counts)
colnames(pb_counts) <- samples

write.csv(
  data.frame(gene = rownames(pb_counts), pb_counts, check.names = FALSE),
  file.path(OBJ_DIR, "expanded_EC_pseudobulk_raw_counts_all_genes.csv"),
  row.names = FALSE
)

pb_meta <- sample_meta %>%
  filter(sample_id %in% colnames(pb_counts)) %>%
  arrange(match(sample_id, colnames(pb_counts)))
pb_meta$condition <- factor(as.character(pb_meta$condition), levels = c("PA", "AC"))
pb_meta$patient_id <- factor(pb_meta$patient_id)
write.csv(pb_meta, file.path(OUT_DIR, "09_pseudobulk_sample_metadata.csv"), row.names = FALSE)

edgeR_available <- requireNamespace("edgeR", quietly = TRUE)
edgeR_success <- FALSE
edgeR_all <- NULL
edgeR_candidates <- NULL

if (edgeR_available) {
  edgeR_success <- tryCatch({
    suppressPackageStartupMessages(library(edgeR))
    design <- model.matrix(~ patient_id + condition, data = pb_meta)
    if (qr(design)$rank < ncol(design)) stop("Paired design matrix is not full rank.")

    y_all <- edgeR::DGEList(counts = pb_counts[, pb_meta$sample_id, drop = FALSE])
    y_all <- edgeR::calcNormFactors(y_all)
    logcpm_all <- edgeR::cpm(y_all, log = TRUE, prior.count = 1)

    pb_target_rows <- lapply(present_targets, function(g) {
      f <- actual_features[[g]]
      data.frame(
        gene = g,
        sample_id = colnames(logcpm_all),
        logCPM = as.numeric(logcpm_all[f, colnames(logcpm_all)]),
        raw_count = as.numeric(pb_counts[f, colnames(logcpm_all)]),
        stringsAsFactors = FALSE
      )
    })
    pb_target_expr <- bind_rows(pb_target_rows) %>%
      left_join(pb_meta %>% select(sample_id, patient_id, condition), by = "sample_id")
    write.csv(pb_target_expr, file.path(OUT_DIR, "10_candidate_pseudobulk_logCPM_and_counts.csv"), row.names = FALSE)

    keep <- edgeR::filterByExpr(y_all, design = design, min.count = 5)
    y <- y_all[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    coef_idx <- which(colnames(design) == "conditionAC")
    if (length(coef_idx) != 1) stop("conditionAC coefficient not found.")
    qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
    tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
    tab$matched_feature <- rownames(tab)
    tab <- tab[, c("matched_feature", setdiff(colnames(tab), "matched_feature"))]
    edgeR_all <<- tab
    write.csv(tab, file.path(OBJ_DIR, "edgeR_paired_AC_vs_PA_all_genes.csv"), row.names = FALSE)

    tested_lookup <- setNames(seq_len(nrow(tab)), toupper(tab$matched_feature))
    candidate_test_rows <- lapply(target_genes, function(g) {
      avail <- gene_availability[gene_availability$gene == g, , drop = FALSE]
      if (!avail$present_in_object) {
        return(data.frame(
          gene = g, role = avail$role, matched_feature = NA_character_,
          edgeR_status = "gene_absent", logFC = NA_real_, logCPM = NA_real_,
          F = NA_real_, PValue = NA_real_, FDR = NA_real_
        ))
      }
      idx2 <- tested_lookup[[toupper(avail$matched_feature)]]
      if (is.null(idx2) || is.na(idx2)) {
        return(data.frame(
          gene = g, role = avail$role, matched_feature = avail$matched_feature,
          edgeR_status = "filtered_low_expression", logFC = NA_real_, logCPM = NA_real_,
          F = NA_real_, PValue = NA_real_, FDR = NA_real_
        ))
      }
      rr <- tab[idx2, , drop = FALSE]
      data.frame(
        gene = g,
        role = avail$role,
        matched_feature = avail$matched_feature,
        edgeR_status = "tested",
        logFC = rr$logFC,
        logCPM = rr$logCPM,
        F = rr$F,
        PValue = rr$PValue,
        FDR = rr$FDR
      )
    })
    edgeR_candidates <<- bind_rows(candidate_test_rows)
    write.csv(edgeR_candidates, file.path(OUT_DIR, "11_candidate_edgeR_paired_AC_vs_PA.csv"), row.names = FALSE)

    pb_stats <- paired_stats(pb_target_expr, "logCPM")
    write.csv(pb_stats, file.path(OUT_DIR, "12_paired_statistics_pseudobulk_logCPM.csv"), row.names = FALSE)

    p_pb <- ggplot(
      pb_target_expr,
      aes(condition, logCPM, group = patient_id, colour = patient_id)
    ) +
      geom_line(linewidth = 0.75, alpha = 0.85) +
      geom_point(size = 2.6) +
      facet_wrap(~ gene, scales = "free_y", ncol = 5) +
      labs(
        x = NULL,
        y = "Pseudobulk logCPM",
        title = "Patient-level pseudobulk: paired PA versus AC"
      )
    save_png(p_pb, "05_candidate_paired_pseudobulk_logCPM.png", 16, 8.5)

    # 每位患者 AC-PA 的 logCPM 差值热图。
    pb_diff <- pb_target_expr %>%
      select(gene, patient_id, condition, logCPM) %>%
      pivot_wider(names_from = condition, values_from = logCPM) %>%
      mutate(AC_minus_PA = AC - PA) %>%
      select(gene, patient_id, AC_minus_PA) %>%
      pivot_wider(names_from = patient_id, values_from = AC_minus_PA)
    pb_diff_mat <- as.matrix(pb_diff[, setdiff(colnames(pb_diff), "gene"), drop = FALSE])
    rownames(pb_diff_mat) <- pb_diff$gene
    grDevices::png(
      file.path(OUT_DIR, "06_patient_AC_minus_PA_logCPM_heatmap.png"),
      width = 2400, height = 2200, res = 300, bg = "white"
    )
    pheatmap::pheatmap(
      pb_diff_mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = NA,
      display_numbers = TRUE,
      number_format = "%.2f",
      main = "Patient-level pseudobulk difference: AC minus PA"
    )
    grDevices::dev.off()

    tested_plot <- edgeR_candidates %>% filter(edgeR_status == "tested")
    if (nrow(tested_plot) > 0) {
      p_edger <- ggplot(tested_plot, aes(x = reorder(gene, logFC), y = logFC)) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(size = -log10(pmax(FDR, 1e-300)))) +
        coord_flip() +
        labs(
          x = NULL,
          y = "edgeR log2 fold change (AC versus PA)",
          size = "-log10 FDR",
          title = "Paired pseudobulk edgeR for nine candidates and SULF1"
        )
      save_png(p_edger, "07_candidate_edgeR_logFC.png", 9, 6.5)
    }

    TRUE
  }, error = function(e) {
    writeLines(
      c("edgeR analysis failed:", conditionMessage(e)),
      file.path(OUT_DIR, "11_edgeR_NOT_RUN.txt")
    )
    FALSE
  })
} else {
  # 没有 edgeR 时仍保存简单 CPM，用于描述性检查。
  lib_size <- Matrix::colSums(pb_counts)
  simple_logcpm <- log2(t(t(pb_counts) / pmax(lib_size, 1)) * 1e6 + 1)
  pb_target_rows <- lapply(present_targets, function(g) {
    f <- actual_features[[g]]
    data.frame(
      gene = g,
      sample_id = colnames(simple_logcpm),
      logCPM = as.numeric(simple_logcpm[f, colnames(simple_logcpm)]),
      raw_count = as.numeric(pb_counts[f, colnames(simple_logcpm)]),
      stringsAsFactors = FALSE
    )
  })
  pb_target_expr <- bind_rows(pb_target_rows) %>%
    left_join(pb_meta %>% select(sample_id, patient_id, condition), by = "sample_id")
  write.csv(pb_target_expr, file.path(OUT_DIR, "10_candidate_pseudobulk_logCPM_and_counts.csv"), row.names = FALSE)
  write.csv(paired_stats(pb_target_expr, "logCPM"), file.path(OUT_DIR, "12_paired_statistics_pseudobulk_logCPM.csv"), row.names = FALSE)
  writeLines(
    "edgeR is not installed. Descriptive pseudobulk logCPM was generated, but formal paired testing was skipped.",
    file.path(OUT_DIR, "11_edgeR_NOT_RUN.txt")
  )
}

# =============================================================================
# E. 汇总：不把“方向一致”夸大为因果或时间顺序证明
# =============================================================================
pb_stats_file <- file.path(OUT_DIR, "12_paired_statistics_pseudobulk_logCPM.csv")
pb_stats_summary <- if (file.exists(pb_stats_file)) read.csv(pb_stats_file, check.names = FALSE) else data.frame()

summary_table <- gene_availability %>%
  left_join(
    avg_stats %>%
      select(
        gene,
        avgexpr_mean_AC_minus_PA = mean_AC_minus_PA,
        avgexpr_n_AC_higher = n_AC_higher,
        avgexpr_n_PA_higher = n_PA_higher,
        avgexpr_direction_pattern = direction_pattern,
        avgexpr_wilcoxon_p = wilcoxon_p,
        avgexpr_paired_Cohen_dz = paired_Cohen_dz
      ),
    by = "gene"
  ) %>%
  left_join(
    pct_stats %>%
      select(
        gene,
        detection_mean_AC_minus_PA = mean_AC_minus_PA,
        detection_n_AC_higher = n_AC_higher,
        detection_n_PA_higher = n_PA_higher,
        detection_direction_pattern = direction_pattern
      ),
    by = "gene"
  )

if (nrow(pb_stats_summary) > 0) {
  summary_table <- summary_table %>%
    left_join(
      pb_stats_summary %>%
        select(
          gene,
          pseudobulk_mean_AC_minus_PA = mean_AC_minus_PA,
          pseudobulk_n_AC_higher = n_AC_higher,
          pseudobulk_n_PA_higher = n_PA_higher,
          pseudobulk_direction_pattern = direction_pattern,
          pseudobulk_wilcoxon_p = wilcoxon_p,
          pseudobulk_paired_Cohen_dz = paired_Cohen_dz
        ),
      by = "gene"
    )
} else {
  summary_table$pseudobulk_mean_AC_minus_PA <- NA_real_
  summary_table$pseudobulk_n_AC_higher <- NA_integer_
  summary_table$pseudobulk_n_PA_higher <- NA_integer_
  summary_table$pseudobulk_direction_pattern <- NA_character_
  summary_table$pseudobulk_wilcoxon_p <- NA_real_
  summary_table$pseudobulk_paired_Cohen_dz <- NA_real_
}

if (!is.null(edgeR_candidates) && nrow(edgeR_candidates) > 0) {
  summary_table <- summary_table %>%
    left_join(
      edgeR_candidates %>%
        select(
          gene,
          edgeR_status,
          edgeR_logFC = logFC,
          edgeR_PValue = PValue,
          edgeR_FDR = FDR
        ),
      by = "gene"
    )
} else {
  summary_table$edgeR_status <- ifelse(summary_table$present_in_object, "not_completed", "gene_absent")
  summary_table$edgeR_logFC <- NA_real_
  summary_table$edgeR_PValue <- NA_real_
  summary_table$edgeR_FDR <- NA_real_
}

summary_table <- summary_table %>%
  mutate(
    statistical_status = case_when(
      !present_in_object ~ "gene_absent",
      edgeR_status == "filtered_low_expression" ~ "descriptive_only_low_expression",
      is.finite(edgeR_FDR) & edgeR_FDR < 0.05 ~ "FDR_below_0.05",
      is.finite(edgeR_PValue) & edgeR_PValue < 0.05 ~ "nominal_P_below_0.05_only",
      edgeR_status == "tested" ~ "tested_not_statistically_significant",
      TRUE ~ "formal_test_not_completed"
    ),
    external_pattern = case_when(
      !present_in_object ~ "not_assessed_gene_absent",
      is.na(pseudobulk_n_AC_higher) | is.na(pseudobulk_n_PA_higher) ~ "pseudobulk_direction_not_available",
      pseudobulk_n_AC_higher == 3 ~ "directionally_consistent_AC_higher_3_of_3",
      pseudobulk_n_PA_higher == 3 ~ "directionally_consistent_PA_higher_3_of_3",
      pseudobulk_n_AC_higher == 2 ~ "partially_consistent_AC_higher_2_of_3",
      pseudobulk_n_PA_higher == 2 ~ "partially_consistent_PA_higher_2_of_3",
      TRUE ~ "mixed_or_tied"
    ),
    interpretation_limit = paste0(
      "PA-AC is a spatial paired contrast; this result does not by itself prove temporal precedence before SULF1 or causal regulation."
    )
  )
write.csv(summary_table, file.path(OUT_DIR, "13_final_candidate_validation_summary.csv"), row.names = FALSE)

# SULF1 与候选的样本级描述性相关性；只有 6 个样本，不作强推断。
if (exists("pb_target_expr") && "SULF1" %in% pb_target_expr$gene) {
  wide_pb <- pb_target_expr %>%
    select(sample_id, gene, logCPM) %>%
    pivot_wider(names_from = gene, values_from = logCPM)
  sul_cor <- lapply(intersect(candidate_genes, colnames(wide_pb)), function(g) {
    rho <- suppressWarnings(cor(wide_pb[[g]], wide_pb[["SULF1"]], method = "spearman", use = "complete.obs"))
    data.frame(
      gene = g,
      spearman_rho_with_SULF1_across_6_samples = rho,
      n_samples = sum(complete.cases(wide_pb[, c(g, "SULF1")])),
      interpretation = "descriptive_only_n_equals_6_not_temporal_or_causal"
    )
  })
  write.csv(bind_rows(sul_cor), file.path(OUT_DIR, "14_candidate_SULF1_sample_correlation_descriptive.csv"), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

n_present_candidates <- sum(gene_availability$present_in_object & gene_availability$role == "early_candidate")
summary_text <- c(
  "GSE159677 nine early-candidate validation completed.",
  paste0("Input expanded EC cells: ", ncol(ec)),
  paste0("Early candidates present in expression matrix: ", n_present_candidates, "/9"),
  paste0("SULF1 anchor present: ", gene_availability$present_in_object[gene_availability$gene == "SULF1"]),
  paste0("Formal paired edgeR completed: ", edgeR_success),
  "Biological replicate: patient, not cell.",
  "PA = plaque-proximal adjacent diseased artery; it is not a truly healthy control.",
  "AC = atherosclerotic core region of the plaque.",
  "AC versus PA is a within-patient spatial contrast, not a longitudinal early-versus-late time series.",
  "A 3/3 paired direction supports cross-patient spatial reproducibility.",
  "It does not prove that the candidate occurs before SULF1, drives SULF1, or causes endothelial transition.",
  "Mapped-state and Tier analyses are supporting sensitivity analyses, not lineage tracing."
)
writeLines(summary_text, file.path(OUT_DIR, "README_9candidate_results.txt"))
cat(paste(summary_text, collapse = "\n"), "\n")
