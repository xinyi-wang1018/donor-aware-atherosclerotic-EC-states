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

# Discovery EC reference already established in the 53-sample study.
REF_RDS <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"
REF_STATE_COL <- "EC_cluster_std"

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

OUT  <- file.path(BASE, "step3B_reference_state_mapping")
TAB  <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ  <- file.path(OUT, "objects")
LOG  <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step3B_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 STEP 3B\n")
cat("Reference EC-state mapping: C7 / C3 / C5 / C1-C2 + C4 off-axis\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

for (fp in c(REF_RDS, MAIN_RDS, EXPANDED_RDS)) {
  if (!file.exists(fp)) stop("Required RDS missing: ", fp)
}

# -------------------------------------------------------------------------
# 0. NON-CIRCULAR MAPPING RULE
# -------------------------------------------------------------------------
# These candidates are NEVER used to:
# - select anchor features
# - calculate reference centroids
# - transfer labels
# - decide mapping confidence
#
# They are examined only AFTER state labels are assigned.
FORBIDDEN_CANDIDATES <- c("ANGPT2", "LOX", "CTHRC1", "SULF1")

# Primary axis state order.
AXIS_STATES <- c(
  "C7_homeostatic",
  "C3_early",
  "C5_transition",
  "C1C2_late"
)

ALL_STATES <- c(
  AXIS_STATES,
  "C4_offaxis",
  "Unassigned"
)

# Conservative pre-specified mapping confidence thresholds.
TRANSFER_SCORE_STRICT <- 0.50
TRANSFER_MARGIN_STRICT <- 0.10
CORRELATION_MIN_STRICT <- 0.10
CORRELATION_MARGIN_STRICT <- 0.03

TRANSFER_SCORE_SENS <- 0.35
TRANSFER_MARGIN_SENS <- 0.05

mapping_rules <- data.frame(
  item = c(
    "Reference RDS",
    "Reference state column",
    "Forbidden candidate genes",
    "Strict transfer score",
    "Strict transfer margin",
    "Strict centroid correlation minimum",
    "Strict centroid correlation margin",
    "Sensitivity transfer score",
    "Sensitivity transfer margin"
  ),
  value = c(
    REF_RDS,
    REF_STATE_COL,
    paste(FORBIDDEN_CANDIDATES, collapse = ";"),
    TRANSFER_SCORE_STRICT,
    TRANSFER_MARGIN_STRICT,
    CORRELATION_MIN_STRICT,
    CORRELATION_MARGIN_STRICT,
    TRANSFER_SCORE_SENS,
    TRANSFER_MARGIN_SENS
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mapping_rules,
  file.path(TAB, "00_mapping_rules_and_forbidden_genes.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 1. Helper functions
# -------------------------------------------------------------------------
join_rna_layers <- function(sobj) {
  DefaultAssay(sobj) <- "RNA"

  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    sobj[["RNA"]] <- tryCatch(
      SeuratObject::JoinLayers(sobj[["RNA"]]),
      error = function(e) sobj[["RNA"]]
    )
  }
  sobj
}

prepare_reference <- function(ref) {
  ref <- join_rna_layers(ref)
  DefaultAssay(ref) <- "RNA"

  if (!REF_STATE_COL %in% colnames(ref@meta.data)) {
    stop("Reference state column missing: ", REF_STATE_COL)
  }

  raw_state <- as.character(ref[[REF_STATE_COL]][,1])

  keep <- raw_state %in% c("C1","C2","C3","C4","C5","C7")
  ref <- subset(ref, cells = colnames(ref)[keep])

  raw_state <- as.character(ref[[REF_STATE_COL]][,1])

  ref$mapping_state <- dplyr::case_when(
    raw_state %in% c("C1","C2") ~ "C1C2_late",
    raw_state == "C3" ~ "C3_early",
    raw_state == "C4" ~ "C4_offaxis",
    raw_state == "C5" ~ "C5_transition",
    raw_state == "C7" ~ "C7_homeostatic",
    TRUE ~ NA_character_
  )

  if (any(is.na(ref$mapping_state))) {
    stop("Unexpected reference state after collapsing.")
  }

  # Cap very large states to reduce abundance-driven label-transfer bias.
  max_per_state <- 800
  keep_cells <- unlist(
    lapply(
      split(colnames(ref), ref$mapping_state),
      function(cc) {
        if (length(cc) <= max_per_state) cc
        else sample(cc, max_per_state)
      }
    ),
    use.names = FALSE
  )

  ref <- subset(ref, cells = keep_cells)

  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(
    ref,
    selection.method = "vst",
    nfeatures = 3500,
    verbose = FALSE
  )

  # Remove forbidden candidates and highly technical/non-specific features.
  features <- VariableFeatures(ref)
  features <- features[
    !features %in% FORBIDDEN_CANDIDATES &
      !grepl("^MT-", features, ignore.case = TRUE) &
      !grepl("^RPL|^RPS", features) &
      !grepl("^HBA|^HBB", features) &
      !grepl("^IG[HKL]", features) &
      features != "MALAT1"
  ]

  list(ref = ref, candidate_features = features)
}

prepare_query <- function(query, tech_keep) {
  query <- join_rna_layers(query)
  DefaultAssay(query) <- "RNA"

  if (!all(c("tech","disease_stage","donor_id") %in%
           colnames(query@meta.data))) {
    stop("Query lacks tech/disease_stage/donor_id metadata.")
  }

  query <- subset(query, subset = tech == tech_keep)

  query$disease_stage <- factor(
    as.character(query$disease_stage),
    levels = c("Normal","Early","Advanced"),
    ordered = TRUE
  )

  query$donor_stage_unit <- paste(
    query$donor_id,
    as.character(query$disease_stage),
    sep = "__"
  )

  query <- NormalizeData(query, verbose = FALSE)
  query
}

score_margin <- function(mat) {
  # mat: cells x states.
  apply(mat, 1, function(x) {
    sx <- sort(x, decreasing = TRUE, na.last = NA)
    if (length(sx) < 2) return(NA_real_)
    sx[1] - sx[2]
  })
}

max_label <- function(mat) {
  apply(mat, 1, function(x) {
    if (all(!is.finite(x))) return(NA_character_)
    colnames(mat)[which.max(x)]
  })
}

map_one_query <- function(ref, candidate_features, query,
                          analysis_set, tech_keep) {

  query <- prepare_query(query, tech_keep)

  shared <- intersect(candidate_features, rownames(query))

  if (length(shared) < 1000) {
    stop(
      "Too few shared non-circular mapping features: ",
      length(shared)
    )
  }

  # Use top shared reference-variable features.
  features <- head(shared, 2500)

  # Explicit check.
  if (any(FORBIDDEN_CANDIDATES %in% features)) {
    stop("Forbidden candidates accidentally entered mapping feature set.")
  }

  cat(
    analysis_set, " / ", tech_keep,
    ": query cells=", ncol(query),
    ", mapping features=", length(features), "\n",
    sep = ""
  )

  # Reference PCA.
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
    reduction.name = "pca",
    verbose = FALSE
  )

  # Query scaling is also needed for the independent centroid method.
  query <- ScaleData(
    query,
    features = features,
    verbose = FALSE
  )

  # ----------------------------
  # Method 1: Seurat label transfer
  # ----------------------------
  anchors <- FindTransferAnchors(
    reference = ref2,
    query = query,
    reference.assay = "RNA",
    query.assay = "RNA",
    normalization.method = "LogNormalize",
    reduction = "pcaproject",
    reference.reduction = "pca",
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

  score_cols <- grep(
    "^prediction.score\\.",
    colnames(pred),
    value = TRUE
  )
  score_cols <- setdiff(score_cols, "prediction.score.max")

  transfer_scores <- as.matrix(pred[, score_cols, drop = FALSE])
  colnames(transfer_scores) <- sub(
    "^prediction.score\\.",
    "",
    colnames(transfer_scores)
  )

  transfer_label <- as.character(pred$predicted.id)
  transfer_max <- as.numeric(pred$prediction.score.max)
  transfer_margin <- score_margin(transfer_scores)

  # ----------------------------
  # Method 2: independent scaled-centroid correlation
  # ----------------------------
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

  state_levels <- c(
    "C7_homeostatic",
    "C3_early",
    "C5_transition",
    "C1C2_late",
    "C4_offaxis"
  )

  centroid <- sapply(state_levels, function(st) {
    cells <- colnames(ref2)[ref2$mapping_state == st]
    Matrix::rowMeans(ref_scaled[, cells, drop = FALSE])
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

  corr_label <- max_label(corr)
  corr_max <- apply(corr, 1, max, na.rm = TRUE)
  corr_margin <- score_margin(corr)

  # Cell order alignment.
  cells <- colnames(query)

  transfer_df <- data.frame(
    cell = pred$cell,
    transfer_label = transfer_label,
    transfer_score = transfer_max,
    transfer_margin = transfer_margin,
    stringsAsFactors = FALSE
  )

  corr_df <- data.frame(
    cell = rownames(corr),
    centroid_label = corr_label,
    centroid_correlation = corr_max,
    centroid_margin = corr_margin,
    stringsAsFactors = FALSE
  )

  audit <- data.frame(
    cell = cells,
    stringsAsFactors = FALSE
  ) %>%
    left_join(transfer_df, by = "cell") %>%
    left_join(corr_df, by = "cell")

  audit <- audit %>%
    mutate(
      methods_agree =
        !is.na(transfer_label) &
        !is.na(centroid_label) &
        transfer_label == centroid_label,

      mapping_call_strict = case_when(
        methods_agree &
          transfer_score >= TRANSFER_SCORE_STRICT &
          transfer_margin >= TRANSFER_MARGIN_STRICT &
          centroid_correlation >= CORRELATION_MIN_STRICT &
          centroid_margin >= CORRELATION_MARGIN_STRICT
          ~ transfer_label,
        TRUE ~ "Unassigned"
      ),

      mapping_call_sensitivity = case_when(
        methods_agree &
          transfer_score >= TRANSFER_SCORE_SENS &
          transfer_margin >= TRANSFER_MARGIN_SENS
          ~ transfer_label,
        TRUE ~ "Unassigned"
      )
    )

  # Attach query metadata.
  audit <- audit %>%
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

  audit$analysis_set <- analysis_set

  # Store predictions in query object.
  audit2 <- audit[match(colnames(query), audit$cell), , drop = FALSE]
  stopifnot(identical(colnames(query), audit2$cell))

  query$mapping_call_strict <- audit2$mapping_call_strict
  query$mapping_call_sensitivity <- audit2$mapping_call_sensitivity
  query$transfer_label <- audit2$transfer_label
  query$transfer_score <- audit2$transfer_score
  query$transfer_margin <- audit2$transfer_margin
  query$centroid_label <- audit2$centroid_label
  query$centroid_correlation <- audit2$centroid_correlation
  query$centroid_margin <- audit2$centroid_margin

  # Add a 2-D query PCA/UMAP only for visualization; it does not affect labels.
  query <- FindVariableFeatures(
    query,
    nfeatures = 2000,
    verbose = FALSE
  )
  q_features <- setdiff(
    intersect(VariableFeatures(query), rownames(query)),
    FORBIDDEN_CANDIDATES
  )
  q_features <- head(q_features, 2000)

  query <- ScaleData(
    query,
    features = q_features,
    verbose = FALSE
  )
  query <- RunPCA(
    query,
    features = q_features,
    npcs = 25,
    reduction.name = "pca.mapview",
    reduction.key = "PCMAP_",
    verbose = FALSE
  )
  query <- RunUMAP(
    query,
    reduction = "pca.mapview",
    dims = 1:20,
    reduction.name = "umap.mapview",
    reduction.key = "UMAPMAP_",
    seed.use = 20260814,
    verbose = FALSE
  )

  list(
    query = query,
    audit = audit,
    features = features,
    reference = ref2
  )
}

make_state_proportions <- function(audit, call_col, analysis_set) {

  x <- audit %>%
    mutate(
      state_call = .data[[call_col]]
    )

  totals <- x %>%
    count(
      donor_id,
      disease_stage,
      donor_stage_unit,
      tech,
      name = "n_total_EC"
    )

  cnt <- x %>%
    count(
      donor_id,
      disease_stage,
      donor_stage_unit,
      tech,
      state_call,
      name = "n_state"
    ) %>%
    complete(
      nesting(
        donor_id,
        disease_stage,
        donor_stage_unit,
        tech
      ),
      state_call = ALL_STATES,
      fill = list(n_state = 0)
    ) %>%
    left_join(
      totals,
      by = c(
        "donor_id",
        "disease_stage",
        "donor_stage_unit",
        "tech"
      )
    ) %>%
    mutate(
      proportion_all_EC = n_state / n_total_EC,
      analysis_set = analysis_set,
      mapping_definition = call_col
    )

  cnt
}

run_state_tests <- function(prop_df, analysis_set) {

  d <- prop_df %>%
    filter(
      analysis_set == .env$analysis_set,
      state_call %in% AXIS_STATES
    )

  out <- lapply(AXIS_STATES, function(st) {
    s <- d %>% filter(state_call == st)

    N <- s %>% filter(as.character(disease_stage) == "Normal")
    E <- s %>% filter(as.character(disease_stage) == "Early")
    A <- s %>% filter(as.character(disease_stage) == "Advanced")

    w_en <- tryCatch(
      wilcox.test(
        E$proportion_all_EC,
        N$proportion_all_EC,
        exact = FALSE
      ),
      error = function(e) NULL
    )

    w_an <- tryCatch(
      wilcox.test(
        A$proportion_all_EC,
        N$proportion_all_EC,
        exact = FALSE
      ),
      error = function(e) NULL
    )

    pair <- inner_join(
      E %>% select(donor_id, E = proportion_all_EC),
      A %>% select(donor_id, A = proportion_all_EC),
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
        contrast = "Early-Normal",
        mean_effect = mean(E$proportion_all_EC) -
          mean(N$proportion_all_EC),
        median_effect = median(E$proportion_all_EC) -
          median(N$proportion_all_EC),
        p_value = ifelse(is.null(w_en), NA, w_en$p.value),
        n_num = nrow(E),
        n_den = nrow(N),
        n_paired = NA_integer_,
        n_pairs_positive = NA_integer_
      ),
      data.frame(
        state = st,
        contrast = "Advanced-Normal",
        mean_effect = mean(A$proportion_all_EC) -
          mean(N$proportion_all_EC),
        median_effect = median(A$proportion_all_EC) -
          median(N$proportion_all_EC),
        p_value = ifelse(is.null(w_an), NA, w_an$p.value),
        n_num = nrow(A),
        n_den = nrow(N),
        n_paired = NA_integer_,
        n_pairs_positive = NA_integer_
      ),
      data.frame(
        state = st,
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
        p_value = ifelse(is.null(w_ea), NA, w_ea$p.value),
        n_num = NA_integer_,
        n_den = NA_integer_,
        n_paired = nrow(pair),
        n_pairs_positive = ifelse(
          nrow(pair) > 0,
          sum(pair$A - pair$E > 0),
          NA_integer_
        )
      )
    )
  }) %>% bind_rows()

  out %>%
    group_by(contrast) %>%
    mutate(state_panel_FDR = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(analysis_set = analysis_set)
}

# Candidate expression AFTER mapping, using donor-stage-state means.
candidate_expression_by_state <- function(mapped_query, analysis_set) {

  DefaultAssay(mapped_query) <- "RNA"
  mapped_query <- NormalizeData(mapped_query, verbose = FALSE)

  data_mat <- GetAssayData(
    mapped_query,
    assay = "RNA",
    layer = "data"
  )

  genes <- intersect(
    FORBIDDEN_CANDIDATES,
    rownames(data_mat)
  )

  md <- mapped_query@meta.data
  md$cell <- rownames(md)

  keep_cells <- rownames(md)[
    md$mapping_call_strict %in%
      c(
        "C7_homeostatic",
        "C3_early",
        "C5_transition",
        "C1C2_late",
        "C4_offaxis"
      )
  ]

  if (length(keep_cells) == 0 || length(genes) == 0) {
    return(data.frame())
  }

  md2 <- md[keep_cells, , drop = FALSE]

  group_id <- paste(
    md2$donor_stage_unit,
    md2$mapping_call_strict,
    sep = "__STATE__"
  )

  groups <- split(rownames(md2), group_id)

  rows <- lapply(names(groups), function(gr) {
    cc <- groups[[gr]]
    first <- md2[cc[1], , drop = FALSE]

    if (length(cc) < 3) return(NULL)

    vals <- Matrix::rowMeans(
      data_mat[genes, cc, drop = FALSE]
    )

    data.frame(
      donor_id = first$donor_id,
      disease_stage = first$disease_stage,
      donor_stage_unit = first$donor_stage_unit,
      mapped_state = first$mapping_call_strict,
      n_cells = length(cc),
      gene = genes,
      mean_log_normalized_expression = as.numeric(vals),
      analysis_set = analysis_set,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(rows)
}

# -------------------------------------------------------------------------
# 2. Prepare reference
# -------------------------------------------------------------------------
ref_raw <- readRDS(REF_RDS)

cat("Discovery reference cells:", ncol(ref_raw), "\n")
cat("Discovery reference genes:", nrow(ref_raw), "\n")
cat(
  "Reference states:",
  paste(sort(unique(as.character(ref_raw[[REF_STATE_COL]][,1]))),
        collapse = ", "),
  "\n\n"
)

prep <- prepare_reference(ref_raw)
ref <- prep$ref
candidate_features <- prep$candidate_features

ref_state_counts <- data.frame(
  state = names(table(ref$mapping_state)),
  n_reference_cells = as.integer(table(ref$mapping_state)),
  stringsAsFactors = FALSE
)

write.csv(
  ref_state_counts,
  file.path(TAB, "01_reference_state_counts_after_balance_cap.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 3. Map MAIN strict FL ECs
# -------------------------------------------------------------------------
main_all <- readRDS(MAIN_RDS)

main_map <- map_one_query(
  ref = ref,
  candidate_features = candidate_features,
  query = main_all,
  analysis_set = "MAIN_A_STRICT_EC",
  tech_keep = "FL"
)

write.csv(
  main_map$audit,
  gzfile(file.path(TAB, "02_MAIN_FL_cell_level_mapping_audit.csv.gz")),
  row.names = FALSE
)

writeLines(
  main_map$features,
  file.path(TAB, "03_MAPPING_FEATURES_MAIN_FL.txt")
)

# Hard proof that candidates were not used.
feature_audit <- data.frame(
  gene = FORBIDDEN_CANDIDATES,
  used_in_mapping = FORBIDDEN_CANDIDATES %in% main_map$features,
  stringsAsFactors = FALSE
)

write.csv(
  feature_audit,
  file.path(TAB, "04_forbidden_candidate_feature_audit.csv"),
  row.names = FALSE
)

if (any(feature_audit$used_in_mapping)) {
  stop("Circularity guard failed: forbidden candidate entered mapping.")
}

# Mapping confidence summary.
main_conf <- main_map$audit %>%
  group_by(disease_stage) %>%
  summarise(
    n_cells = n(),
    pct_methods_agree = mean(methods_agree) * 100,
    pct_strict_assigned =
      mean(mapping_call_strict != "Unassigned") * 100,
    pct_sensitivity_assigned =
      mean(mapping_call_sensitivity != "Unassigned") * 100,
    median_transfer_score = median(transfer_score, na.rm = TRUE),
    median_transfer_margin = median(transfer_margin, na.rm = TRUE),
    median_centroid_correlation =
      median(centroid_correlation, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  main_conf,
  file.path(TAB, "05_MAIN_FL_mapping_confidence_by_stage.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 4. State proportions and donor-level stage tests
# -------------------------------------------------------------------------
main_prop_strict <- make_state_proportions(
  main_map$audit,
  "mapping_call_strict",
  "MAIN_A_STRICT_EC"
)

main_prop_senscall <- make_state_proportions(
  main_map$audit,
  "mapping_call_sensitivity",
  "MAIN_A_STRICT_EC_sensitivity_mapping_threshold"
)

write.csv(
  main_prop_strict,
  file.path(TAB, "06_MAIN_FL_donor_stage_state_proportions_STRICT.csv"),
  row.names = FALSE
)

write.csv(
  main_prop_senscall,
  file.path(TAB, "07_MAIN_FL_donor_stage_state_proportions_SENSITIVITY_THRESHOLD.csv"),
  row.names = FALSE
)

main_state_tests <- run_state_tests(
  main_prop_strict,
  "MAIN_A_STRICT_EC"
)

write.csv(
  main_state_tests,
  file.path(TAB, "08_MAIN_FL_state_proportion_stage_tests.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 5. Expanded A+B FL sensitivity mapping
# -------------------------------------------------------------------------
expanded_all <- readRDS(EXPANDED_RDS)

expanded_map <- map_one_query(
  ref = ref,
  candidate_features = candidate_features,
  query = expanded_all,
  analysis_set = "EXPANDED_AB_EC",
  tech_keep = "FL"
)

write.csv(
  expanded_map$audit,
  gzfile(file.path(TAB, "09_EXPANDED_FL_cell_level_mapping_audit.csv.gz")),
  row.names = FALSE
)

expanded_prop <- make_state_proportions(
  expanded_map$audit,
  "mapping_call_strict",
  "EXPANDED_AB_EC"
)

write.csv(
  expanded_prop,
  file.path(TAB, "10_EXPANDED_FL_donor_stage_state_proportions_STRICT.csv"),
  row.names = FALSE
)

expanded_state_tests <- run_state_tests(
  expanded_prop,
  "EXPANDED_AB_EC"
)

write.csv(
  expanded_state_tests,
  file.path(TAB, "11_EXPANDED_FL_state_proportion_stage_tests.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 6. MAIN vs EXPANDED directional concordance
# -------------------------------------------------------------------------
concordance <- main_state_tests %>%
  select(
    state,
    contrast,
    main_mean_effect = mean_effect,
    main_median_effect = median_effect,
    main_P = p_value,
    main_FDR = state_panel_FDR
  ) %>%
  left_join(
    expanded_state_tests %>%
      select(
        state,
        contrast,
        expanded_mean_effect = mean_effect,
        expanded_median_effect = median_effect,
        expanded_P = p_value,
        expanded_FDR = state_panel_FDR
      ),
    by = c("state","contrast")
  ) %>%
  mutate(
    same_direction =
      sign(main_mean_effect) == sign(expanded_mean_effect)
  )

write.csv(
  concordance,
  file.path(TAB, "12_state_proportion_MAIN_vs_EXPANDED_concordance.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 7. Non-circular check: excluded candidate genes AFTER mapping
# -------------------------------------------------------------------------
main_candidate_by_state <- candidate_expression_by_state(
  main_map$query,
  "MAIN_A_STRICT_EC"
)

expanded_candidate_by_state <- candidate_expression_by_state(
  expanded_map$query,
  "EXPANDED_AB_EC"
)

candidate_by_state <- bind_rows(
  main_candidate_by_state,
  expanded_candidate_by_state
)

write.csv(
  candidate_by_state,
  file.path(TAB, "13_EXCLUDED_candidates_expression_by_mapped_state.csv"),
  row.names = FALSE
)

candidate_state_summary <- candidate_by_state %>%
  group_by(
    analysis_set,
    mapped_state,
    gene
  ) %>%
  summarise(
    n_donor_stage_state_units = n(),
    mean_expression =
      mean(mean_log_normalized_expression),
    median_expression =
      median(mean_log_normalized_expression),
    .groups = "drop"
  )

write.csv(
  candidate_state_summary,
  file.path(TAB, "14_EXCLUDED_candidates_mapped_state_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8. 3DE strict EC mapping - DESCRIPTIVE ONLY
# -------------------------------------------------------------------------
de3_map <- map_one_query(
  ref = ref,
  candidate_features = candidate_features,
  query = main_all,
  analysis_set = "MAIN_A_STRICT_EC_3DE",
  tech_keep = "3DE"
)

de3_prop <- make_state_proportions(
  de3_map$audit,
  "mapping_call_strict",
  "MAIN_A_STRICT_EC_3DE"
)

write.csv(
  de3_prop,
  file.path(TAB, "15_3DE_state_proportions_DESCRIPTIVE_ONLY.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 9. Evidence classification
# -------------------------------------------------------------------------
get_effect <- function(tbl, st, contrast) {
  x <- tbl %>%
    filter(state == st, .data$contrast == contrast)
  if (nrow(x) == 0) return(NA_real_)
  x$mean_effect[1]
}

evidence <- data.frame(
  biological_question = c(
    "C7 decreases by Early",
    "C7 decreases by Advanced",
    "C3 enriched in Early",
    "C1C2 enriched in Advanced vs Normal",
    "C1C2 enriched in Advanced vs matched Early"
  ),
  state = c(
    "C7_homeostatic",
    "C7_homeostatic",
    "C3_early",
    "C1C2_late",
    "C1C2_late"
  ),
  contrast = c(
    "Early-Normal",
    "Advanced-Normal",
    "Early-Normal",
    "Advanced-Normal",
    "Advanced-Early_paired"
  ),
  expected_sign = c(-1,-1,1,1,1),
  stringsAsFactors = FALSE
)

evidence$main_effect <- mapply(
  function(st, ct) get_effect(main_state_tests, st, ct),
  evidence$state,
  evidence$contrast
)

evidence$expanded_effect <- mapply(
  function(st, ct) get_effect(expanded_state_tests, st, ct),
  evidence$state,
  evidence$contrast
)

evidence$main_expected_direction <- with(
  evidence,
  sign(main_effect) == expected_sign
)

evidence$expanded_expected_direction <- with(
  evidence,
  sign(expanded_effect) == expected_sign
)

evidence$evidence_label <- dplyr::case_when(
  evidence$main_expected_direction &
    evidence$expanded_expected_direction
    ~ "direction_reproduced_in_MAIN_and_EXPANDED",

  evidence$main_expected_direction
    ~ "direction_supported_in_MAIN_only",

  TRUE
    ~ "expected_stage_enrichment_not_reproduced"
)

write.csv(
  evidence,
  file.path(TAB, "16_primary_state_composition_evidence_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 10. Figures - PURE WHITE BACKGROUND
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

# UMAP of query mapping labels.
p1 <- DimPlot(
  main_map$query,
  reduction = "umap.mapview",
  group.by = "mapping_call_strict",
  raster = TRUE,
  pt.size = 0.45,
  shuffle = TRUE
) +
  ggtitle("GSE218166 FL strict EC: reference-state mapping") +
  theme_classic(base_size = 10)

save_white(
  p1,
  "01_MAIN_FL_reference_state_mapping_UMAP.png",
  11,
  8
)

# UMAP split by pathology stage - visualization only.
p2 <- DimPlot(
  main_map$query,
  reduction = "umap.mapview",
  group.by = "mapping_call_strict",
  split.by = "disease_stage",
  raster = TRUE,
  pt.size = 0.35,
  ncol = 3
) +
  ggtitle("Mapped EC states across pathology stages") +
  theme_classic(base_size = 9)

save_white(
  p2,
  "02_MAIN_FL_mapped_states_split_by_pathology_stage.png",
  16,
  6
)

# Donor-stage state proportions.
plot_prop <- main_prop_strict %>%
  filter(state_call %in% ALL_STATES)

p3 <- ggplot(
  plot_prop,
  aes(
    x = disease_stage,
    y = proportion_all_EC
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.2) +
  facet_wrap(
    ~state_call,
    scales = "free_y",
    ncol = 3
  ) +
  theme_classic(base_size = 10) +
  labs(
    title = "Donor-level EC-state proportions by pathology stage",
    subtitle = "Denominator = all final strict ECs in each donor-stage unit",
    x = NULL,
    y = "Proportion of ECs"
  )

save_white(
  p3,
  "03_MAIN_FL_donor_state_proportions_by_stage.png",
  13,
  8
)

# Stacked mean composition - descriptive.
mean_comp <- main_prop_strict %>%
  group_by(disease_stage, state_call) %>%
  summarise(
    mean_proportion = mean(proportion_all_EC),
    .groups = "drop"
  )

p4 <- ggplot(
  mean_comp,
  aes(
    x = disease_stage,
    y = mean_proportion,
    fill = state_call
  )
) +
  geom_col() +
  theme_classic(base_size = 10) +
  labs(
    title = "Mean EC-state composition across pathology stages",
    subtitle = "Descriptive; inferential tests are donor-level",
    x = NULL,
    y = "Mean proportion",
    fill = "Mapped state"
  )

save_white(
  p4,
  "04_MAIN_FL_mean_state_composition_stacked.png",
  9,
  7
)

# Paired Early-Advanced changes in state proportions.
paired_donors <- intersect(
  unique(
    main_prop_strict$donor_id[
      as.character(main_prop_strict$disease_stage) == "Early"
    ]
  ),
  unique(
    main_prop_strict$donor_id[
      as.character(main_prop_strict$disease_stage) == "Advanced"
    ]
  )
)

paired_prop <- main_prop_strict %>%
  filter(
    donor_id %in% paired_donors,
    as.character(disease_stage) %in% c("Early","Advanced"),
    state_call %in% AXIS_STATES
  )

p5 <- ggplot(
  paired_prop,
  aes(
    x = disease_stage,
    y = proportion_all_EC,
    group = donor_id
  )
) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~state_call,
    scales = "free_y",
    ncol = 2
  ) +
  theme_classic(base_size = 10) +
  labs(
    title = "Matched-donor Early to Advanced EC-state changes",
    subtitle = paste(
      "Paired donors:",
      paste(paired_donors, collapse = ", ")
    ),
    x = NULL,
    y = "Proportion of strict ECs"
  )

save_white(
  p5,
  "05_MAIN_FL_paired_Early_Advanced_state_proportions.png",
  11,
  8
)

# Mapping confidence.
conf_long <- main_map$audit %>%
  select(
    disease_stage,
    transfer_score,
    centroid_correlation
  ) %>%
  pivot_longer(
    cols = c(transfer_score, centroid_correlation),
    names_to = "metric",
    values_to = "value"
  )

p6 <- ggplot(
  conf_long,
  aes(
    x = disease_stage,
    y = value
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 0.8) +
  facet_wrap(~metric, scales = "free_y") +
  theme_classic(base_size = 10) +
  labs(
    title = "Reference-mapping confidence by pathology stage",
    x = NULL,
    y = NULL
  )

save_white(
  p6,
  "06_MAIN_FL_mapping_confidence_by_stage.png",
  10,
  6
)

# Excluded candidates AFTER non-circular state mapping.
if (nrow(main_candidate_by_state) > 0) {
  p7 <- ggplot(
    main_candidate_by_state,
    aes(
      x = mapped_state,
      y = mean_log_normalized_expression
    )
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 1.7) +
    facet_wrap(
      ~gene,
      scales = "free_y",
      ncol = 2
    ) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1)
    ) +
    labs(
      title = "Excluded candidate genes AFTER reference-state mapping",
      subtitle = "ANGPT2/LOX/CTHRC1/SULF1 were not used to assign states",
      x = NULL,
      y = "Donor-stage-state mean log-normalized expression"
    )

  save_white(
    p7,
    "07_NONCIRCULAR_candidate_expression_by_mapped_state.png",
    13,
    8
  )
}

# -------------------------------------------------------------------------
# 11. Save mapped objects on server
# -------------------------------------------------------------------------
saveRDS(
  main_map$query,
  file.path(
    OBJ,
    "GSE218166_STEP3B_MAIN_FL_strict_EC_mapped.rds"
  ),
  compress = FALSE
)

saveRDS(
  expanded_map$query,
  file.path(
    OBJ,
    "GSE218166_STEP3B_EXPANDED_FL_EC_mapped.rds"
  ),
  compress = FALSE
)

saveRDS(
  de3_map$query,
  file.path(
    OBJ,
    "GSE218166_STEP3B_3DE_strict_EC_mapped_DESCRIPTIVE.rds"
  ),
  compress = FALSE
)

# -------------------------------------------------------------------------
# 12. Human-readable summary
# -------------------------------------------------------------------------
strict_state_counts <- main_map$audit %>%
  count(disease_stage, mapping_call_strict, name = "n_cells") %>%
  group_by(disease_stage) %>%
  mutate(
    fraction = n_cells / sum(n_cells)
  ) %>%
  ungroup()

strict_assigned <- main_map$audit %>%
  summarise(
    n_total = n(),
    n_assigned = sum(mapping_call_strict != "Unassigned"),
    pct_assigned =
      100 * mean(mapping_call_strict != "Unassigned")
  )

evidence_lines <- evidence %>%
  mutate(
    line = paste0(
      biological_question,
      " | MAIN effect=",
      sprintf("%.3f", main_effect),
      " | EXPANDED effect=",
      sprintf("%.3f", expanded_effect),
      " | ",
      evidence_label
    )
  ) %>%
  pull(line)

summary_lines <- c(
  "============================================================",
  "GSE218166 STEP3B REFERENCE STATE MAPPING",
  "============================================================",
  paste0("Run time: ", format(Sys.time())),
  "",
  "REFERENCE",
  paste0("RDS: ", REF_RDS),
  paste0("State column: ", REF_STATE_COL),
  "Reference C1+C2 were collapsed as C1C2_late.",
  "C4 was retained as an off-axis competitor to avoid forcing every EC onto the main trajectory.",
  "",
  "NON-CIRCULARITY GUARD",
  paste0(
    "Excluded from ALL mapping features: ",
    paste(FORBIDDEN_CANDIDATES, collapse = ", ")
  ),
  paste0(
    "Forbidden genes accidentally used: ",
    ifelse(any(feature_audit$used_in_mapping), "YES", "NO")
  ),
  "",
  "MAPPING METHOD",
  "A state is called in the strict analysis only when Seurat label transfer and an independent scaled-centroid correlation agree,",
  "and both confidence/margin thresholds are passed.",
  "Otherwise the cell is Unassigned rather than being forced into a trajectory state.",
  "",
  paste0("MAIN FL strict EC total: ", strict_assigned$n_total),
  paste0(
    "Strictly assigned: ",
    strict_assigned$n_assigned,
    " (",
    sprintf("%.1f", strict_assigned$pct_assigned),
    "%)"
  ),
  "",
  "PRIMARY BIOLOGICAL QUESTION",
  "Does pathology progression change the COMPOSITION of EC transcriptional states rather than making all ECs change synchronously?",
  "",
  "DIRECTIONAL STATE-COMPOSITION SUMMARY",
  evidence_lines,
  "",
  "STATISTICAL UNIT",
  "State proportions are calculated per donor x pathology-stage unit.",
  "Normal-vs-Early and Normal-vs-Advanced are donor-level unpaired tests.",
  "Early-vs-Advanced uses only donors observed in both stages.",
  "No cell-level P values are used for stage-composition inference.",
  "",
  "3DE",
  "Mapped for descriptive technology sensitivity only; no equal-strength Early-stage inference because the final strict-EC Early group has one donor.",
  "",
  "NEXT",
  "Upload GSE218166_step3B_REVIEW_SMALL.tar.gz.",
  "Then decide whether the discovery trajectory is externally supported as state-composition redistribution, partially supported, or not reproduced."
)

writeLines(
  summary_lines,
  file.path(OUT, "Step3B_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\nCompleted:", format(Sys.time()), "\n")
