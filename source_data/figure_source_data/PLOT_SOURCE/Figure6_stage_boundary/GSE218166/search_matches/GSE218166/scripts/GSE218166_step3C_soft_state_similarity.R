#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260814)

BASE <- "/root/GSE218166"

REF_RDS <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"

MAIN_RDS <- file.path(
  BASE,
  "step2E_FINAL_AUTHOR_ALIGNED_EC_FREEZE",
  "objects",
  "GSE218166_FINAL_MAIN_A_STRICT_EC.rds"
)

EXPANDED_RDS <- file.path(
  BASE,
  "step2E_FINAL_AUTHOR_ALIGNED_EC_FREEZE",
  "objects",
  "GSE218166_FINAL_EXPANDED_AB_EC.rds"
)

OUT  <- file.path(BASE, "step3C_soft_state_similarity")
TAB  <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ  <- file.path(OUT, "objects")
LOG  <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step3C_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 STEP 3C - SOFT REFERENCE-STATE SIMILARITY\n")
cat("All final ECs contribute continuous state-similarity scores.\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

for (fp in c(REF_RDS, MAIN_RDS, EXPANDED_RDS)) {
  if (!file.exists(fp)) stop("Required file missing: ", fp)
}

FORBIDDEN <- c("ANGPT2", "LOX", "CTHRC1", "SULF1")
STATE_COL_CANDIDATES <- c("EC_cluster_std", "EC_cluster")

STATE_LEVELS <- c(
  "C7_homeostatic",
  "C3_early",
  "C5_transition",
  "C1C2_late",
  "C4_offaxis"
)

# -------------------------------------------------------------------------
# 1. Reference-state handling
# -------------------------------------------------------------------------
canonicalize_state <- function(v) {
  x <- trimws(as.character(v))
  numeric_like <- grepl("^[0-9]+$", x)
  x[numeric_like] <- paste0("C", x[numeric_like])
  x <- sub("^c([0-9]+)$", "C\\1", x, ignore.case = TRUE)
  x
}

detect_state_col <- function(ref) {
  hit <- STATE_COL_CANDIDATES[
    STATE_COL_CANDIDATES %in% colnames(ref@meta.data)
  ]
  if (length(hit) == 0) {
    stop("No EC_cluster_std or EC_cluster column in reference.")
  }
  hit[1]
}

join_layers_safe <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj[["RNA"]] <- tryCatch(
      SeuratObject::JoinLayers(obj[["RNA"]]),
      error = function(e) obj[["RNA"]]
    )
  }
  obj
}

prepare_reference <- function(ref) {

  ref <- join_layers_safe(ref)
  DefaultAssay(ref) <- "RNA"

  state_col <- detect_state_col(ref)
  raw <- canonicalize_state(ref[[state_col]][,1])

  required <- c("C1","C2","C3","C4","C5","C7")
  missing <- setdiff(required, unique(raw))
  if (length(missing) > 0) {
    stop("Reference is missing states: ", paste(missing, collapse = ", "))
  }

  keep <- raw %in% required
  ref <- subset(ref, cells = colnames(ref)[keep])
  raw <- canonicalize_state(ref[[state_col]][,1])

  ref$mapping_state <- case_when(
    raw %in% c("C1","C2") ~ "C1C2_late",
    raw == "C3" ~ "C3_early",
    raw == "C4" ~ "C4_offaxis",
    raw == "C5" ~ "C5_transition",
    raw == "C7" ~ "C7_homeostatic",
    TRUE ~ NA_character_
  )

  # Balance only very large classes; preserve all smaller classes.
  max_per_state <- 800L
  cells_keep <- unlist(
    lapply(
      split(colnames(ref), ref$mapping_state),
      function(cc) {
        if (length(cc) <= max_per_state) cc
        else sample(cc, max_per_state)
      }
    ),
    use.names = FALSE
  )

  ref <- subset(ref, cells = cells_keep)
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(
    ref,
    selection.method = "vst",
    nfeatures = 3500,
    verbose = FALSE
  )

  features <- VariableFeatures(ref)
  features <- features[
    !features %in% FORBIDDEN &
      !grepl("^MT-", features, ignore.case = TRUE) &
      !grepl("^RPL|^RPS", features) &
      !grepl("^HBA|^HBB", features) &
      !grepl("^IG[HKL]", features) &
      features != "MALAT1"
  ]

  list(
    ref = ref,
    state_col = state_col,
    features = features
  )
}

prepare_query <- function(obj, tech_keep = "FL") {

  obj <- join_layers_safe(obj)
  DefaultAssay(obj) <- "RNA"

  obj <- subset(obj, subset = tech == tech_keep)

  obj$disease_stage <- factor(
    as.character(obj$disease_stage),
    levels = c("Normal","Early","Advanced"),
    ordered = TRUE
  )

  obj$donor_stage_unit <- paste(
    obj$donor_id,
    as.character(obj$disease_stage),
    sep = "__"
  )

  obj <- NormalizeData(obj, verbose = FALSE)
  obj
}

# -------------------------------------------------------------------------
# 2. Soft mapping: full transfer probabilities + full centroid similarities
# -------------------------------------------------------------------------
soft_map <- function(ref, candidate_features, query,
                     analysis_set, tech_keep = "FL") {

  query <- prepare_query(query, tech_keep)

  shared <- intersect(candidate_features, rownames(query))
  features <- head(shared, 2500)

  if (length(features) < 1000) {
    stop("Too few shared mapping genes: ", length(features))
  }

  if (any(FORBIDDEN %in% features)) {
    stop("Circularity guard failed: forbidden gene in mapping features.")
  }

  cat(
    analysis_set, " / ", tech_keep,
    ": cells=", ncol(query),
    ", mapping features=", length(features), "\n",
    sep = ""
  )

  ref2 <- ref

  ref2 <- ScaleData(
    ref2,
    features = features,
    verbose = FALSE
  )

  ref2 <- RunPCA(
    ref2,
    features = features,
    npcs = 30,
    reduction.name = "pca.softmap",
    reduction.key = "PCSOFT_",
    verbose = FALSE
  )

  query <- ScaleData(
    query,
    features = features,
    verbose = FALSE
  )

  anchors <- FindTransferAnchors(
    reference = ref2,
    query = query,
    reference.assay = "RNA",
    query.assay = "RNA",
    normalization.method = "LogNormalize",
    reduction = "pcaproject",
    reference.reduction = "pca.softmap",
    features = features,
    dims = 1:30,
    verbose = FALSE
  )

  pred <- TransferData(
    anchorset = anchors,
    refdata = ref2$mapping_state,
    dims = 1:30,
    verbose = FALSE
  )

  pred$cell <- rownames(pred)

  transfer_cols <- grep(
    "^prediction.score\\.",
    colnames(pred),
    value = TRUE
  )
  transfer_cols <- setdiff(
    transfer_cols,
    "prediction.score.max"
  )

  transfer <- pred[, c("cell", transfer_cols), drop = FALSE]
  colnames(transfer) <- c(
    "cell",
    sub("^prediction.score\\.", "", transfer_cols)
  )

  transfer_long <- transfer %>%
    pivot_longer(
      cols = -cell,
      names_to = "state",
      values_to = "transfer_probability"
    )

  # Independent centroid correlations.
  ref_scaled <- as.matrix(
    GetAssayData(
      ref2,
      assay = "RNA",
      layer = "scale.data"
    )[features, , drop = FALSE]
  )

  query_scaled <- as.matrix(
    GetAssayData(
      query,
      assay = "RNA",
      layer = "scale.data"
    )[features, , drop = FALSE]
  )

  centroid <- sapply(STATE_LEVELS, function(st) {
    cc <- colnames(ref2)[ref2$mapping_state == st]
    Matrix::rowMeans(ref_scaled[, cc, drop = FALSE])
  })

  centroid <- as.matrix(centroid)
  rownames(centroid) <- features

  corr <- suppressWarnings(
    cor(
      query_scaled,
      centroid,
      use = "pairwise.complete.obs",
      method = "pearson"
    )
  )

  corr <- as.matrix(corr)

  corr_df <- as.data.frame(corr)
  corr_df$cell <- rownames(corr_df)

  corr_long <- corr_df %>%
    pivot_longer(
      cols = -cell,
      names_to = "state",
      values_to = "centroid_correlation"
    )

  # Join both independent continuous similarity measures.
  long <- transfer_long %>%
    full_join(
      corr_long,
      by = c("cell","state")
    ) %>%
    left_join(
      query@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        select(
          cell,
          donor_id,
          disease_stage,
          donor_stage_unit,
          tech
        ),
      by = "cell"
    )

  long$analysis_set <- analysis_set

  list(
    query = query,
    long = long,
    features = features
  )
}

# -------------------------------------------------------------------------
# 3. Donor-stage aggregation
# -------------------------------------------------------------------------
aggregate_soft <- function(long_df) {

  long_df %>%
    group_by(
      analysis_set,
      donor_id,
      disease_stage,
      donor_stage_unit,
      tech,
      state
    ) %>%
    summarise(
      n_cells = n(),
      mean_transfer_probability =
        mean(transfer_probability, na.rm = TRUE),
      median_transfer_probability =
        median(transfer_probability, na.rm = TRUE),
      mean_centroid_correlation =
        mean(centroid_correlation, na.rm = TRUE),
      median_centroid_correlation =
        median(centroid_correlation, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------------------------------------------------------------------------
# 4. Donor-level stage tests for continuous state similarity
# -------------------------------------------------------------------------
test_one_metric <- function(agg, metric_col, analysis_set) {

  states <- STATE_LEVELS

  out <- lapply(states, function(st) {

    d <- agg %>%
      filter(
        .data$analysis_set == .env$analysis_set,
        state == st
      )

    N <- d %>% filter(as.character(disease_stage) == "Normal")
    E <- d %>% filter(as.character(disease_stage) == "Early")
    A <- d %>% filter(as.character(disease_stage) == "Advanced")

    Nval <- N[[metric_col]]
    Eval <- E[[metric_col]]
    Aval <- A[[metric_col]]

    w_en <- tryCatch(
      wilcox.test(Eval, Nval, exact = FALSE),
      error = function(e) NULL
    )

    w_an <- tryCatch(
      wilcox.test(Aval, Nval, exact = FALSE),
      error = function(e) NULL
    )

    pair <- inner_join(
      E %>% select(donor_id, E = all_of(metric_col)),
      A %>% select(donor_id, A = all_of(metric_col)),
      by = "donor_id"
    )

    w_ea <- if (nrow(pair) >= 3) {
      tryCatch(
        wilcox.test(
          pair$A,
          pair$E,
          paired = TRUE,
          exact = FALSE
        ),
        error = function(e) NULL
      )
    } else NULL

    bind_rows(
      data.frame(
        state = st,
        metric = metric_col,
        contrast = "Early-Normal",
        mean_effect = mean(Eval) - mean(Nval),
        median_effect = median(Eval) - median(Nval),
        p_value = ifelse(is.null(w_en), NA, w_en$p.value),
        n_num = length(Eval),
        n_den = length(Nval),
        n_paired = NA_integer_
      ),
      data.frame(
        state = st,
        metric = metric_col,
        contrast = "Advanced-Normal",
        mean_effect = mean(Aval) - mean(Nval),
        median_effect = median(Aval) - median(Nval),
        p_value = ifelse(is.null(w_an), NA, w_an$p.value),
        n_num = length(Aval),
        n_den = length(Nval),
        n_paired = NA_integer_
      ),
      data.frame(
        state = st,
        metric = metric_col,
        contrast = "Advanced-Early_paired",
        mean_effect = ifelse(
          nrow(pair) > 0,
          mean(pair$A - pair$E),
          NA_real_
        ),
        median_effect = ifelse(
          nrow(pair) > 0,
          median(pair$A - pair$E),
          NA_real_
        ),
        p_value = ifelse(
          is.null(w_ea),
          NA,
          w_ea$p.value
        ),
        n_num = NA_integer_,
        n_den = NA_integer_,
        n_paired = nrow(pair)
      )
    )
  }) %>% bind_rows()

  out %>%
    group_by(metric, contrast) %>%
    mutate(state_FDR = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(analysis_set = analysis_set)
}

run_tests <- function(agg, analysis_set) {
  bind_rows(
    test_one_metric(
      agg,
      "mean_transfer_probability",
      analysis_set
    ),
    test_one_metric(
      agg,
      "mean_centroid_correlation",
      analysis_set
    )
  )
}

# -------------------------------------------------------------------------
# 5. Run MAIN + EXPANDED FL
# -------------------------------------------------------------------------
ref_raw <- readRDS(REF_RDS)
prep <- prepare_reference(ref_raw)

cat("Reference state column:", prep$state_col, "\n")
cat("Reference cells after cap:", ncol(prep$ref), "\n")
print(table(prep$ref$mapping_state))
cat("\n")

main_all <- readRDS(MAIN_RDS)
expanded_all <- readRDS(EXPANDED_RDS)

main_map <- soft_map(
  prep$ref,
  prep$features,
  main_all,
  "MAIN_A_STRICT_EC",
  "FL"
)

expanded_map <- soft_map(
  prep$ref,
  prep$features,
  expanded_all,
  "EXPANDED_AB_EC",
  "FL"
)

writeLines(
  main_map$features,
  file.path(TAB, "01_SOFT_MAPPING_FEATURES.txt")
)

forbidden_audit <- data.frame(
  gene = FORBIDDEN,
  used_in_mapping = FORBIDDEN %in% main_map$features,
  stringsAsFactors = FALSE
)

write.csv(
  forbidden_audit,
  file.path(TAB, "02_FORBIDDEN_CANDIDATE_AUDIT.csv"),
  row.names = FALSE
)

write.csv(
  main_map$long,
  gzfile(file.path(TAB, "03_MAIN_FL_cell_soft_state_scores.csv.gz")),
  row.names = FALSE
)

write.csv(
  expanded_map$long,
  gzfile(file.path(TAB, "04_EXPANDED_FL_cell_soft_state_scores.csv.gz")),
  row.names = FALSE
)

main_agg <- aggregate_soft(main_map$long)
expanded_agg <- aggregate_soft(expanded_map$long)

write.csv(
  main_agg,
  file.path(TAB, "05_MAIN_FL_donor_stage_soft_state_scores.csv"),
  row.names = FALSE
)

write.csv(
  expanded_agg,
  file.path(TAB, "06_EXPANDED_FL_donor_stage_soft_state_scores.csv"),
  row.names = FALSE
)

main_tests <- run_tests(
  main_agg,
  "MAIN_A_STRICT_EC"
)

expanded_tests <- run_tests(
  expanded_agg,
  "EXPANDED_AB_EC"
)

write.csv(
  main_tests,
  file.path(TAB, "07_MAIN_FL_soft_state_stage_tests.csv"),
  row.names = FALSE
)

write.csv(
  expanded_tests,
  file.path(TAB, "08_EXPANDED_FL_soft_state_stage_tests.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 6. Four-way concordance:
# MAIN / EXPANDED x Transfer / Centroid
# -------------------------------------------------------------------------
main_wide <- main_tests %>%
  select(
    state,
    metric,
    contrast,
    main_effect = mean_effect,
    main_P = p_value,
    main_FDR = state_FDR
  )

expanded_wide <- expanded_tests %>%
  select(
    state,
    metric,
    contrast,
    expanded_effect = mean_effect,
    expanded_P = p_value,
    expanded_FDR = state_FDR
  )

concordance <- main_wide %>%
  left_join(
    expanded_wide,
    by = c("state","metric","contrast")
  ) %>%
  mutate(
    MAIN_EXPANDED_same_direction =
      sign(main_effect) == sign(expanded_effect)
  )

write.csv(
  concordance,
  file.path(TAB, "09_soft_state_MAIN_vs_EXPANDED_concordance.csv"),
  row.names = FALSE
)

# Biological expectation table.
expectations <- data.frame(
  question = c(
    "C7 similarity decreases in Early",
    "C7 similarity decreases in Advanced",
    "C3 similarity increases in Early",
    "C5 transition similarity increases in Early",
    "C1C2 similarity increases in Advanced vs Normal",
    "C1C2 similarity increases in Advanced vs matched Early"
  ),
  state = c(
    "C7_homeostatic",
    "C7_homeostatic",
    "C3_early",
    "C5_transition",
    "C1C2_late",
    "C1C2_late"
  ),
  contrast = c(
    "Early-Normal",
    "Advanced-Normal",
    "Early-Normal",
    "Early-Normal",
    "Advanced-Normal",
    "Advanced-Early_paired"
  ),
  expected_sign = c(-1,-1,1,1,1,1),
  stringsAsFactors = FALSE
)

evidence_rows <- lapply(seq_len(nrow(expectations)), function(i) {

  q <- expectations[i, ]

  x <- concordance %>%
    filter(
      state == q$state,
      contrast == q$contrast
    )

  if (nrow(x) == 0) {
    return(data.frame(
      question = q$question,
      state = q$state,
      contrast = q$contrast,
      expected_sign = q$expected_sign,
      transfer_MAIN = NA,
      transfer_EXPANDED = NA,
      centroid_MAIN = NA,
      centroid_EXPANDED = NA,
      all_four_expected_direction = FALSE,
      evidence_label = "not_evaluable"
    ))
  }

  trow <- x %>% filter(metric == "mean_transfer_probability")
  crow <- x %>% filter(metric == "mean_centroid_correlation")

  vals <- c(
    trow$main_effect[1],
    trow$expanded_effect[1],
    crow$main_effect[1],
    crow$expanded_effect[1]
  )

  ok <- all(sign(vals) == q$expected_sign)

  data.frame(
    question = q$question,
    state = q$state,
    contrast = q$contrast,
    expected_sign = q$expected_sign,
    transfer_MAIN = vals[1],
    transfer_EXPANDED = vals[2],
    centroid_MAIN = vals[3],
    centroid_EXPANDED = vals[4],
    all_four_expected_direction = ok,
    evidence_label = ifelse(
      ok,
      "direction_consistent_across_MAIN_EXPANDED_and_two_methods",
      "not_consistent_across_all_four_checks"
    ),
    stringsAsFactors = FALSE
  )
})

evidence <- bind_rows(evidence_rows)

write.csv(
  evidence,
  file.path(TAB, "10_primary_soft_state_evidence_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 7. 3DE descriptive only
# -------------------------------------------------------------------------
de3_map <- soft_map(
  prep$ref,
  prep$features,
  main_all,
  "MAIN_A_STRICT_EC_3DE",
  "3DE"
)

de3_agg <- aggregate_soft(de3_map$long)

write.csv(
  de3_agg,
  file.path(TAB, "11_3DE_donor_stage_soft_state_scores_DESCRIPTIVE_ONLY.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8. Figures - white background
# -------------------------------------------------------------------------
save_white <- function(p, filename, width = 10, height = 7) {
  ggsave(
    file.path(PLOT, filename),
    p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

# Transfer probabilities.
p1 <- ggplot(
  main_agg,
  aes(
    x = disease_stage,
    y = mean_transfer_probability
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.1) +
  facet_wrap(~state, scales = "free_y", ncol = 3) +
  theme_classic(base_size = 10) +
  labs(
    title = "Soft reference-state similarity: Seurat transfer probabilities",
    subtitle = "Each point = one donor × pathology-stage EC pseudobulk unit",
    x = NULL,
    y = "Mean transfer probability"
  )

save_white(
  p1,
  "01_MAIN_FL_soft_transfer_scores_by_stage.png",
  13,
  8
)

# Centroid correlations.
p2 <- ggplot(
  main_agg,
  aes(
    x = disease_stage,
    y = mean_centroid_correlation
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.1) +
  facet_wrap(~state, scales = "free_y", ncol = 3) +
  theme_classic(base_size = 10) +
  labs(
    title = "Soft reference-state similarity: independent centroid correlations",
    subtitle = "All final strict ECs contribute; no hard assignment threshold",
    x = NULL,
    y = "Mean centroid correlation"
  )

save_white(
  p2,
  "02_MAIN_FL_soft_centroid_scores_by_stage.png",
  13,
  8
)

# Paired Early/Advanced transfer scores.
paired_donors <- intersect(
  unique(
    main_agg$donor_id[
      as.character(main_agg$disease_stage) == "Early"
    ]
  ),
  unique(
    main_agg$donor_id[
      as.character(main_agg$disease_stage) == "Advanced"
    ]
  )
)

paired <- main_agg %>%
  filter(
    donor_id %in% paired_donors,
    as.character(disease_stage) %in% c("Early","Advanced")
  )

p3 <- ggplot(
  paired,
  aes(
    x = disease_stage,
    y = mean_transfer_probability,
    group = donor_id
  )
) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2.1) +
  facet_wrap(~state, scales = "free_y", ncol = 3) +
  theme_classic(base_size = 10) +
  labs(
    title = "Matched-donor Early to Advanced soft state similarity",
    x = NULL,
    y = "Mean transfer probability"
  )

save_white(
  p3,
  "03_MAIN_FL_paired_Early_Advanced_soft_transfer.png",
  13,
  8
)

# Heatmap of donor-stage transfer probabilities.
heat <- main_agg %>%
  select(
    donor_stage_unit,
    disease_stage,
    state,
    mean_transfer_probability
  )

sample_order <- main_agg %>%
  distinct(
    donor_stage_unit,
    donor_id,
    disease_stage
  ) %>%
  arrange(disease_stage, donor_id) %>%
  pull(donor_stage_unit)

heat$donor_stage_unit <- factor(
  heat$donor_stage_unit,
  levels = sample_order
)

p4 <- ggplot(
  heat,
  aes(
    x = donor_stage_unit,
    y = state,
    fill = mean_transfer_probability
  )
) +
  geom_tile() +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "Donor-stage soft state similarity matrix",
    fill = "Transfer\nprobability"
  )

save_white(
  p4,
  "04_MAIN_FL_donor_stage_soft_state_heatmap.png",
  13,
  6.5
)

# Effect concordance plot.
effect_plot <- concordance %>%
  filter(
    contrast %in%
      c(
        "Early-Normal",
        "Advanced-Normal",
        "Advanced-Early_paired"
      )
  ) %>%
  select(
    state,
    metric,
    contrast,
    main_effect,
    expanded_effect
  ) %>%
  pivot_longer(
    cols = c(main_effect, expanded_effect),
    names_to = "analysis_set",
    values_to = "effect"
  )

p5 <- ggplot(
  effect_plot,
  aes(
    x = contrast,
    y = effect,
    shape = analysis_set
  )
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 2.4
  ) +
  facet_grid(metric ~ state, scales = "free_y") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1)
  ) +
  labs(
    title = "Soft-state effects: MAIN vs EXPANDED EC definitions",
    x = NULL,
    y = "Stage effect"
  )

save_white(
  p5,
  "05_soft_state_effect_concordance.png",
  16,
  9
)

# -------------------------------------------------------------------------
# 9. Save compact products
# -------------------------------------------------------------------------
saveRDS(
  list(
    main_donor_stage_scores = main_agg,
    expanded_donor_stage_scores = expanded_agg,
    main_tests = main_tests,
    expanded_tests = expanded_tests,
    evidence = evidence,
    mapping_features = main_map$features
  ),
  file.path(
    OBJ,
    "GSE218166_step3C_soft_state_products.rds"
  ),
  compress = "gzip"
)

# -------------------------------------------------------------------------
# 10. Summary
# -------------------------------------------------------------------------
evidence_lines <- evidence %>%
  mutate(
    line = paste0(
      question,
      " | transfer MAIN=",
      sprintf("%.4f", transfer_MAIN),
      " | transfer EXP=",
      sprintf("%.4f", transfer_EXPANDED),
      " | centroid MAIN=",
      sprintf("%.4f", centroid_MAIN),
      " | centroid EXP=",
      sprintf("%.4f", centroid_EXPANDED),
      " | ",
      evidence_label
    )
  ) %>%
  pull(line)

summary_lines <- c(
  "============================================================",
  "GSE218166 STEP3C SOFT REFERENCE-STATE SIMILARITY",
  "============================================================",
  paste0("Run time: ", format(Sys.time())),
  "",
  "WHY THIS ANALYSIS",
  "Step3B strict hard mapping assigned only 16.1% of MAIN FL ECs.",
  "Step3C therefore uses continuous state-similarity scores for ALL final ECs,",
  "without lowering hard-label thresholds or forcing cells into a discrete state.",
  "",
  "NON-CIRCULARITY",
  paste0(
    "Excluded from mapping: ",
    paste(FORBIDDEN, collapse = ", ")
  ),
  paste0(
    "Forbidden features used: ",
    ifelse(any(forbidden_audit$used_in_mapping), "YES", "NO")
  ),
  "",
  "TWO INDEPENDENT SOFT MEASURES",
  "1. Seurat full prediction probability for each reference state.",
  "2. Independent query-to-reference centroid Pearson correlation for each state.",
  "",
  "BIOLOGICAL EVIDENCE",
  evidence_lines,
  "",
  "INTERPRETATION RULE",
  "A stage-composition signal is considered robust only when its expected direction is consistent",
  "in MAIN strict EC, EXPANDED A+B EC, Seurat transfer probability, and centroid correlation.",
  "Statistical significance is supportive but is not used to override direction discordance.",
  "",
  "3DE",
  "Computed descriptively only because Early has one donor after strict EC freezing.",
  "",
  "NEXT",
  "Upload GSE218166_step3C_REVIEW_SMALL.tar.gz.",
  "This will finalize whether GSE218166 supports the discovery trajectory fully, partially, or only at the state-marker level."
)

writeLines(
  summary_lines,
  file.path(OUT, "Step3C_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\nCompleted:", format(Sys.time()), "\n")
