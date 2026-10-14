options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 20 * 1024^3)

# ============================================================
# STRICT-SOURCE TWO-COHORT EC SLINGSHOT
#
# Strict sources:
#   GSE155512
#   GSE253904 study
#
# Main questions:
# 1. Does C7 -> C3 -> C5 -> C1/C2 remain after strict
#    source restriction?
# 2. Does the result remain in strict arterial-only EC?
# 3. Does C4 remain off the main arterial continuum?
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(SingleCellExperiment)
  library(slingshot)
  library(S4Vectors)
  library(ggplot2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

CORE_IN <- "/root/AS_EC_slingshot/AS_EC_arterial_dysfunction_slingshot_scored.rds"

ARTERIAL_IN <- paste0(
  "/root/AS_EC_arterial_only_C3_sensitivity/01_objects/",
  "AS_EC_strict_arterial_only_reclustered_slingshot.rds"
)

OUT <- "/root/AS_strict_source_2cohort_slingshot"

OBJDIR <- file.path(OUT, "objects")
TABDIR <- file.path(OUT, "tables")
PLOTDIR <- file.path(OUT, "plots")
LOGDIR <- file.path(OUT, "logs")

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper: convert metadata safely to atomic character
# ------------------------------------------------------------

atomic_chr <- function(x) {

  if (is.factor(x)) {
    return(as.character(x))
  }

  if (is.list(x)) {

    return(
      vapply(
        x,
        function(z) {
          if (length(z) == 0 || all(is.na(z))) {
            return(NA_character_)
          }

          paste(
            as.character(z),
            collapse = ";"
          )
        },
        FUN.VALUE = character(1)
      )
    )
  }

  as.character(x)
}

# ------------------------------------------------------------
# Detect metadata column
# ------------------------------------------------------------

find_col <- function(obj, candidates) {

  meta_cols <- colnames(obj@meta.data)

  hit <- candidates[
    candidates %in% meta_cols
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

# ------------------------------------------------------------
# Canonicalize dataset
# ------------------------------------------------------------

canonical_dataset <- function(x) {

  x <- atomic_chr(x)

  out <- x

  out[
    grepl(
      "GSE155512",
      x,
      ignore.case = TRUE
    )
  ] <- "GSE155512"

  # GSE253904 is the parent SuperSeries.
  # If child-series labels occur, treat them as the same study.
  out[
    grepl(
      "GSE253902|GSE253903|GSE253904",
      x,
      ignore.case = TRUE
    )
  ] <- "GSE253904"

  out
}

# ------------------------------------------------------------
# Canonicalize EC state
# ------------------------------------------------------------

canonical_state <- function(obj) {

  state_col <- find_col(
    obj,
    c(
      "EC_cluster",
      "seurat_clusters"
    )
  )

  if (is.na(state_col)) {
    stop(
      "Cannot find EC_cluster or seurat_clusters."
    )
  }

  z <- atomic_chr(
    obj@meta.data[[state_col]]
  )

  z <- trimws(z)

  state <- rep(
    NA_character_,
    length(z)
  )

  already_c <- grepl(
    "^C[0-9]+$",
    z
  )

  state[already_c] <- z[already_c]

  numeric_only <- grepl(
    "^[0-9]+$",
    z
  )

  state[numeric_only] <- paste0(
    "C",
    z[numeric_only]
  )

  state
}

# ------------------------------------------------------------
# Check input
# ------------------------------------------------------------

if (!file.exists(CORE_IN)) {
  stop(
    "CORE input missing: ",
    CORE_IN
  )
}

cat(
  "\n============================================================\n"
)

cat(
  "STRICT-SOURCE TWO-COHORT SLINGSHOT\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Core input:\n",
  CORE_IN,
  "\n\n"
)

# ------------------------------------------------------------
# Read object
# ------------------------------------------------------------

core0 <- readRDS(CORE_IN)

ds_col <- find_col(
  core0,
  c(
    "dataset_id",
    "dataset",
    "study",
    "orig.ident"
  )
)

if (is.na(ds_col)) {
  stop(
    "No dataset metadata column found."
  )
}

core0$strict_dataset <- canonical_dataset(
  core0@meta.data[[ds_col]]
)

core0$strict_state <- canonical_state(
  core0
)

# ------------------------------------------------------------
# Initial audit
# ------------------------------------------------------------

audit_before <- as.data.frame(
  table(
    dataset = core0$strict_dataset,
    state = core0$strict_state,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit_before,
  file.path(
    TABDIR,
    "01_CORE_before_source_filter_dataset_state.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Strict source filter
# ------------------------------------------------------------

STRICT_DATASETS <- c(
  "GSE155512",
  "GSE253904"
)

keep <- (
  core0$strict_dataset %in%
    STRICT_DATASETS
) &
(
  core0$strict_state %in%
    c(
      "C1",
      "C2",
      "C3",
      "C4",
      "C5",
      "C7"
    )
)

strict_cells <- Cells(core0)[keep]

if (length(strict_cells) == 0) {
  stop(
    "No strict-source cells retained."
  )
}

core <- subset(
  core0,
  cells = strict_cells
)

rm(core0)
gc()

cat(
  "Strict-source core cells:",
  ncol(core),
  "\n"
)

cat(
  "Datasets:\n"
)

print(
  table(core$strict_dataset)
)

cat(
  "\nStates:\n"
)

print(
  table(core$strict_state)
)

# ------------------------------------------------------------
# Required state check
# ------------------------------------------------------------

required_main_states <- c(
  "C7",
  "C3",
  "C5",
  "C1",
  "C2"
)

missing_states <- setdiff(
  required_main_states,
  unique(core$strict_state)
)

if (length(missing_states) > 0) {

  stop(
    "Required main-axis states missing after source filtering: ",
    paste(
      missing_states,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Save count tables
# ------------------------------------------------------------

core_counts <- as.data.frame(
  table(
    dataset = core$strict_dataset,
    state = core$strict_state
  ),
  stringsAsFactors = FALSE
)

write.csv(
  core_counts,
  file.path(
    TABDIR,
    "02_CORE_STRICT_SOURCE_dataset_state_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Sample-level composition for audit only
# ------------------------------------------------------------

sample_col <- find_col(
  core,
  c(
    "donor_id",
    "patient_id",
    "patient",
    "donor",
    "Individual",
    "sample_id",
    "orig.ident"
  )
)

if (!is.na(sample_col)) {

  core$strict_bio_unit <- atomic_chr(
    core@meta.data[[sample_col]]
  )

  sample_tab <- unique(
    data.frame(
      dataset = core$strict_dataset,
      biological_unit = core$strict_bio_unit,
      state = core$strict_state,
      stringsAsFactors = FALSE
    )
  )

  write.csv(
    sample_tab,
    file.path(
      TABDIR,
      "03_CORE_STRICT_SOURCE_biounit_state_presence.csv"
    ),
    row.names = FALSE
  )
}

# ============================================================
# Recompute local EC geometry
# ============================================================

reprocess_ec <- function(
    obj,
    tag
) {

  cat(
    "\n============================================================\n"
  )

  cat(
    "REPROCESS:",
    tag,
    "\n"
  )

  cat(
    "============================================================\n"
  )

  if (!"RNA" %in% Assays(obj)) {
    stop(
      tag,
      ": RNA assay not found."
    )
  }

  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = min(
      3000,
      nrow(obj)
    ),
    verbose = FALSE
  )

  vf <- VariableFeatures(obj)

  if (length(vf) < 200) {
    stop(
      tag,
      ": too few variable features."
    )
  }

  obj <- ScaleData(
    obj,
    features = vf,
    verbose = FALSE
  )

  max_pc <- min(
    30,
    length(vf) - 1,
    ncol(obj) - 1
  )

  if (max_pc < 10) {
    stop(
      tag,
      ": too few PCs can be computed."
    )
  }

  obj <- RunPCA(
    obj,
    features = vf,
    npcs = max_pc,
    reduction.name = "pca_strict",
    reduction.key = "PCSTRICT_",
    verbose = FALSE
  )

  dims_use <- seq_len(
    min(
      20,
      ncol(
        Embeddings(
          obj,
          "pca_strict"
        )
      )
    )
  )

  # ----------------------------------------------------------
  # Harmony: try modern API, then older API
  # ----------------------------------------------------------

  harmony_error_1 <- NULL

  obj_h <- tryCatch(
    {
      harmony::RunHarmony(
        object = obj,
        group.by.vars = "strict_dataset",
        reduction = "pca_strict",
        dims.use = dims_use,
        reduction.save = "harmony_strict",
        project.dim = FALSE,
        verbose = FALSE
      )
    },
    error = function(e) {

      harmony_error_1 <<- conditionMessage(e)

      NULL
    }
  )

  if (is.null(obj_h)) {

    cat(
      "Modern Harmony call failed; trying legacy API.\n"
    )

    cat(
      "First Harmony error:",
      harmony_error_1,
      "\n"
    )

    obj_h <- tryCatch(
      {
        harmony::RunHarmony(
          object = obj,
          group.by.vars = "strict_dataset",
          reduction.use = "pca_strict",
          dims.use = dims_use,
          reduction.save = "harmony_strict",
          project.dim = FALSE,
          verbose = FALSE
        )
      },
      error = function(e) {

        stop(
          "Harmony failed under both APIs. Legacy error: ",
          conditionMessage(e)
        )
      }
    )
  }

  obj <- obj_h

  if (!"harmony_strict" %in% Reductions(obj)) {
    stop(
      tag,
      ": harmony_strict reduction missing."
    )
  }

  hdim <- ncol(
    Embeddings(
      obj,
      "harmony_strict"
    )
  )

  h_dims <- seq_len(
    min(
      20,
      hdim
    )
  )

  nn <- min(
    30,
    max(
      5,
      ncol(obj) - 1
    )
  )

  obj <- RunUMAP(
    obj,
    reduction = "harmony_strict",
    dims = h_dims,
    n.neighbors = nn,
    reduction.name = "umap_strict",
    reduction.key = "UMAPSTRICT_",
    verbose = FALSE
  )

  attr(
    obj,
    "strict_harmony_dims"
  ) <- h_dims

  obj
}

# ------------------------------------------------------------
# Reprocess full six-state strict-source CORE object
# ------------------------------------------------------------

core <- reprocess_ec(
  core,
  "CORE_STRICT_SOURCE"
)

saveRDS(
  core,
  file.path(
    OBJDIR,
    "CORE_STRICT_SOURCE_reprocessed.rds"
  )
)

# ============================================================
# Slingshot runner
# ============================================================

run_sling <- function(
    obj,
    states,
    tag,
    end_clusters = NULL
) {

  use_cells <- Cells(obj)[
    obj$strict_state %in%
      states
  ]

  x <- subset(
    obj,
    cells = use_cells
  )

  emb_all <- Embeddings(
    x,
    "harmony_strict"
  )

  dims_use <- seq_len(
    min(
      20,
      ncol(emb_all)
    )
  )

  emb <- emb_all[
    ,
    dims_use,
    drop = FALSE
  ]

  st <- atomic_chr(
    x$strict_state
  )

  if (!"C7" %in% st) {
    stop(
      tag,
      ": C7 missing."
    )
  }

  dummy <- matrix(
    0,
    nrow = 1,
    ncol = nrow(emb)
  )

  rownames(dummy) <- "dummy"
  colnames(dummy) <- rownames(emb)

  sce <- SingleCellExperiment(
    assays = list(
      dummy = dummy
    )
  )

  colData(sce)$state <- factor(
    st,
    levels = states
  )

  reducedDim(
    sce,
    "STRICT_HARMONY"
  ) <- emb

  if (is.null(end_clusters)) {

    sce <- slingshot(
      sce,
      clusterLabels = "state",
      reducedDim = "STRICT_HARMONY",
      start.clus = "C7"
    )

  } else {

    sce <- slingshot(
      sce,
      clusterLabels = "state",
      reducedDim = "STRICT_HARMONY",
      start.clus = "C7",
      end.clus = end_clusters
    )
  }

  lins <- slingLineages(sce)

  pt <- slingPseudotime(
    sce,
    na = TRUE
  )

  wt <- slingCurveWeights(sce)

  # ----------------------------------------------------------
  # Lineage table
  # ----------------------------------------------------------

  lineage_df <- data.frame(
    lineage = names(lins),
    path = vapply(
      lins,
      function(z) {
        paste(
          z,
          collapse = " -> "
        )
      },
      FUN.VALUE = character(1)
    ),
    stringsAsFactors = FALSE
  )

  if (
    is.null(names(lins)) ||
    any(names(lins) == "")
  ) {

    lineage_df$lineage <- paste0(
      "Lineage",
      seq_along(lins)
    )
  }

  write.csv(
    lineage_df,
    file.path(
      TABDIR,
      paste0(
        tag,
        "_01_lineages.csv"
      )
    ),
    row.names = FALSE
  )

  cat(
    "\n",
    tag,
    " lineages:\n",
    sep = ""
  )

  print(lineage_df)

  # ----------------------------------------------------------
  # Cell pseudotime
  # ----------------------------------------------------------

  cell_df <- data.frame(
    cell = rownames(emb),
    state = st,
    dataset = atomic_chr(
      x$strict_dataset
    ),
    stringsAsFactors = FALSE
  )

  if (
    "strict_bio_unit" %in%
    colnames(x@meta.data)
  ) {

    cell_df$biological_unit <- atomic_chr(
      x$strict_bio_unit
    )
  }

  for (j in seq_len(ncol(pt))) {

    nm <- colnames(pt)[j]

    if (
      is.null(nm) ||
      is.na(nm) ||
      nm == ""
    ) {
      nm <- paste0(
        "Lineage",
        j
      )
    }

    cell_df[[
        paste0(
          "pseudotime_",
          nm
        )
      ]
    ] <- pt[, j]

    cell_df[[
        paste0(
          "weight_",
          nm
        )
      ]
    ] <- wt[, j]
  }

  write.csv(
    cell_df,
    file.path(
      TABDIR,
      paste0(
        tag,
        "_02_cell_pseudotime.csv"
      )
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # State pseudotime summaries
  # ----------------------------------------------------------

  summary_rows <- list()
  k <- 1

  for (j in seq_len(ncol(pt))) {

    ln <- colnames(pt)[j]

    if (
      is.null(ln) ||
      is.na(ln) ||
      ln == ""
    ) {
      ln <- paste0(
        "Lineage",
        j
      )
    }

    for (ss in states) {

      vals <- pt[
        st == ss,
        j
      ]

      vals <- vals[
        is.finite(vals)
      ]

      summary_rows[[k]] <- data.frame(
        lineage = ln,
        state = ss,
        n_finite = length(vals),
        median_pseudotime = if (
          length(vals) > 0
        ) {
          median(vals)
        } else {
          NA_real_
        },
        mean_pseudotime = if (
          length(vals) > 0
        ) {
          mean(vals)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )

      k <- k + 1
    }
  }

  state_summary <- do.call(
    rbind,
    summary_rows
  )

  write.csv(
    state_summary,
    file.path(
      TABDIR,
      paste0(
        tag,
        "_03_state_pseudotime_summary.csv"
      )
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # Add pseudotime to Seurat subset
  # ----------------------------------------------------------

  for (j in seq_len(ncol(pt))) {

    ln <- colnames(pt)[j]

    if (
      is.null(ln) ||
      is.na(ln) ||
      ln == ""
    ) {
      ln <- paste0(
        "Lineage",
        j
      )
    }

    x[[
        paste0(
          "strict_",
          tag,
          "_",
          ln,
          "_pseudotime"
        )
      ]
    ] <- pt[, j]
  }

  saveRDS(
    x,
    file.path(
      OBJDIR,
      paste0(
        tag,
        "_slingshot.rds"
      )
    )
  )

  # ----------------------------------------------------------
  # Pseudotime distribution plot
  # ----------------------------------------------------------

  if (ncol(pt) >= 1) {

    first_pt <- pt[, 1]

    plot_df <- data.frame(
      state = factor(
        st,
        levels = states
      ),
      pseudotime = first_pt
    )

    p <- ggplot(
      plot_df,
      aes(
        x = state,
        y = pseudotime
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.62
      ) +
      geom_jitter(
        width = 0.15,
        alpha = 0.25,
        size = 0.45
      ) +
      theme_bw(
        base_size = 12
      ) +
      labs(
        title = paste0(
          tag,
          " — first inferred lineage"
        ),
        x = "EC state",
        y = "Slingshot pseudotime"
      ) +
      theme(
        panel.grid.minor = element_blank()
      )

    ggsave(
      filename = file.path(
        PLOTDIR,
        paste0(
          tag,
          "_pseudotime_by_state.png"
        )
      ),
      plot = p,
      width = 7,
      height = 5,
      dpi = 300,
      bg = "white"
    )
  }

  list(
    object = x,
    sce = sce,
    lineages = lineage_df,
    state_summary = state_summary
  )
}

# ============================================================
# ANALYSIS 1
# Main five-state trajectory, fully data-driven except root C7
# C4 is deliberately excluded from the main arterial continuum.
# ============================================================

MAIN_STATES <- c(
  "C7",
  "C3",
  "C5",
  "C1",
  "C2"
)

main_free <- run_sling(
  obj = core,
  states = MAIN_STATES,
  tag = "CORE_MAIN_FREE",
  end_clusters = NULL
)

# ============================================================
# ANALYSIS 2
# Endpoint-constrained sensitivity only.
# This is NOT the primary inference.
# ============================================================

main_constrained <- run_sling(
  obj = core,
  states = MAIN_STATES,
  tag = "CORE_MAIN_CONSTRAINED_C1_C2",
  end_clusters = c(
    "C1",
    "C2"
  )
)

# ============================================================
# ANALYSIS 3
# Six-state data-driven analysis including C4
# Used only to audit whether C4 stays off-axis.
# ============================================================

SIX_STATES <- c(
  "C7",
  "C3",
  "C5",
  "C1",
  "C2",
  "C4"
)

six_free <- run_sling(
  obj = core,
  states = SIX_STATES,
  tag = "CORE_SIX_STATE_FREE",
  end_clusters = NULL
)

# ============================================================
# Strict-source UMAP plot
# ============================================================

um <- Embeddings(
  core,
  "umap_strict"
)

um_df <- data.frame(
  UMAP1 = um[, 1],
  UMAP2 = um[, 2],
  state = factor(
    core$strict_state,
    levels = c(
      "C7",
      "C3",
      "C5",
      "C1",
      "C2",
      "C4"
    )
  ),
  dataset = core$strict_dataset,
  stringsAsFactors = FALSE
)

p_umap <- ggplot(
  um_df,
  aes(
    UMAP1,
    UMAP2,
    color = state
  )
) +
  geom_point(
    size = 0.65,
    alpha = 0.70
  ) +
  theme_bw(
    base_size = 12
  ) +
  labs(
    title = "Strict-source EC",
    subtitle = "GSE155512 + GSE253904",
    color = "State"
  ) +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  file.path(
    PLOTDIR,
    "CORE_STRICT_SOURCE_UMAP_states.png"
  ),
  p_umap,
  width = 7,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

# ============================================================
# ARTERIAL-ONLY independent sensitivity
# ============================================================

arterial_completed <- FALSE

if (file.exists(ARTERIAL_IN)) {

  cat(
    "\n============================================================\n"
  )

  cat(
    "ARTERIAL-ONLY STRICT-SOURCE SENSITIVITY\n"
  )

  cat(
    "============================================================\n"
  )

  art0 <- readRDS(
    ARTERIAL_IN
  )

  art_ds_col <- find_col(
    art0,
    c(
      "dataset_id",
      "dataset",
      "study",
      "orig.ident"
    )
  )

  if (is.na(art_ds_col)) {
    stop(
      "Arterial object has no dataset column."
    )
  }

  art0$strict_dataset <- canonical_dataset(
    art0@meta.data[[art_ds_col]]
  )

  art0$strict_state <- canonical_state(
    art0
  )

  art_keep <- (
    art0$strict_dataset %in%
      STRICT_DATASETS
  ) &
  (
    art0$strict_state %in%
      MAIN_STATES
  )

  art_cells <- Cells(art0)[
    art_keep
  ]

  if (length(art_cells) < 50) {

    warning(
      "Too few arterial-only strict-source cells; skipping."
    )

  } else {

    art <- subset(
      art0,
      cells = art_cells
    )

    rm(art0)
    gc()

    cat(
      "Arterial strict-source cells:",
      ncol(art),
      "\n"
    )

    print(
      table(
        art$strict_dataset,
        art$strict_state
      )
    )

    art_counts <- as.data.frame(
      table(
        dataset = art$strict_dataset,
        state = art$strict_state
      ),
      stringsAsFactors = FALSE
    )

    write.csv(
      art_counts,
      file.path(
        TABDIR,
        "04_ARTERIAL_STRICT_SOURCE_dataset_state_counts.csv"
      ),
      row.names = FALSE
    )

    art <- reprocess_ec(
      art,
      "ARTERIAL_STRICT_SOURCE"
    )

    saveRDS(
      art,
      file.path(
        OBJDIR,
        "ARTERIAL_STRICT_SOURCE_reprocessed.rds"
      )
    )

    arterial_free <- run_sling(
      obj = art,
      states = MAIN_STATES,
      tag = "ARTERIAL_MAIN_FREE",
      end_clusters = NULL
    )

    arterial_constrained <- run_sling(
      obj = art,
      states = MAIN_STATES,
      tag = "ARTERIAL_MAIN_CONSTRAINED_C1_C2",
      end_clusters = c(
        "C1",
        "C2"
      )
    )

    arterial_completed <- TRUE
  }

} else {

  cat(
    "\nArterial input not found; arterial rerun skipped:\n",
    ARTERIAL_IN,
    "\n"
  )
}

# ============================================================
# Automated interpretation helpers
# ============================================================

path_contains_in_order <- function(
    path_string,
    required
) {

  p <- trimws(
    unlist(
      strsplit(
        path_string,
        "->",
        fixed = TRUE
      )
    )
  )

  idx <- match(
    required,
    p
  )

  if (any(is.na(idx))) {
    return(FALSE)
  }

  all(
    diff(idx) > 0
  )
}

check_expected_main <- function(lineage_df) {

  if (
    is.null(lineage_df) ||
    nrow(lineage_df) == 0
  ) {
    return(FALSE)
  }

  checks <- vapply(
    lineage_df$path,
    function(p) {

      prefix_ok <- path_contains_in_order(
        p,
        c(
          "C7",
          "C3",
          "C5"
        )
      )

      terminal_ok <- grepl(
        "C1$|C2$",
        trimws(p)
      )

      prefix_ok && terminal_ok
    },
    FUN.VALUE = logical(1)
  )

  any(checks)
}

core_expected_supported <- check_expected_main(
  main_free$lineages
)

arterial_expected_supported <- NA

if (arterial_completed) {

  arterial_expected_supported <- check_expected_main(
    arterial_free$lineages
  )
}

# ============================================================
# Final checkpoints
# ============================================================

checkpoints <- data.frame(
  checkpoint = c(
    "strict_source_core_cells",
    "strict_source_core_datasets",
    "C7_present",
    "C3_present",
    "C5_present",
    "C1_present",
    "C2_present",
    "core_free_expected_prefix_supported",
    "arterial_analysis_completed",
    "arterial_free_expected_prefix_supported"
  ),
  value = c(
    as.character(
      ncol(core)
    ),
    paste(
      sort(
        unique(
          core$strict_dataset
        )
      ),
      collapse = ";"
    ),
    as.character(
      "C7" %in%
        core$strict_state
    ),
    as.character(
      "C3" %in%
        core$strict_state
    ),
    as.character(
      "C5" %in%
        core$strict_state
    ),
    as.character(
      "C1" %in%
        core$strict_state
    ),
    as.character(
      "C2" %in%
        core$strict_state
    ),
    as.character(
      core_expected_supported
    ),
    as.character(
      arterial_completed
    ),
    as.character(
      arterial_expected_supported
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  checkpoints,
  file.path(
    TABDIR,
    "20_FINAL_CHECKPOINTS.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Final report
# ============================================================

report_file <- file.path(
  OUT,
  "00_FINAL_STRICT_SOURCE_SLINGSHOT_REPORT.txt"
)

sink(
  report_file
)

cat(
  "============================================================\n"
)

cat(
  "STRICT-SOURCE TWO-COHORT EC SLINGSHOT REPORT\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Strict sources:\n"
)

cat(
  "GSE155512 + GSE253904 study\n\n"
)

cat(
  "Core strict-source cells:",
  ncol(core),
  "\n\n"
)

cat(
  "Dataset composition:\n"
)

print(
  table(
    core$strict_dataset
  )
)

cat(
  "\nState composition:\n"
)

print(
  table(
    core$strict_state
  )
)

cat(
  "\n============================================================\n"
)

cat(
  "PRIMARY DATA-DRIVEN FIVE-STATE SLINGSHOT\n"
)

cat(
  "============================================================\n"
)

print(
  main_free$lineages
)

cat(
  "\nState pseudotime summaries:\n"
)

print(
  main_free$state_summary
)

cat(
  "\nExpected C7 -> C3 -> C5 prefix with C1/C2 terminal supported:",
  core_expected_supported,
  "\n"
)

cat(
  "\n============================================================\n"
)

cat(
  "C1/C2 ENDPOINT-CONSTRAINED SENSITIVITY\n"
)

cat(
  "============================================================\n"
)

print(
  main_constrained$lineages
)

cat(
  "\n============================================================\n"
)

cat(
  "SIX-STATE DATA-DRIVEN AUDIT INCLUDING C4\n"
)

cat(
  "============================================================\n"
)

print(
  six_free$lineages
)

cat(
  "\n============================================================\n"
)

cat(
  "ARTERIAL-ONLY STRICT-SOURCE RESULT\n"
)

cat(
  "============================================================\n"
)

cat(
  "Completed:",
  arterial_completed,
  "\n"
)

if (arterial_completed) {

  print(
    arterial_free$lineages
  )

  cat(
    "\nExpected C7 -> C3 -> C5 prefix with C1/C2 terminal supported:",
    arterial_expected_supported,
    "\n"
  )

  cat(
    "\nArterial state pseudotime summaries:\n"
  )

  print(
    arterial_free$state_summary
  )
}

cat(
  "\n============================================================\n"
)

cat(
  "INTERPRETATION RULE\n"
)

cat(
  "============================================================\n"
)

cat(
  paste0(
    "Primary support requires the unforced five-state Slingshot ",
    "to retain C7 -> C3 -> C5 before reaching C1 and/or C2. ",
    "The endpoint-constrained model is sensitivity analysis only. ",
    "The six-state model is used to inspect whether C4 remains ",
    "off the main arterial continuum. The arterial-only analysis ",
    "tests whether the inferred organization persists after ",
    "stricter arterial restriction.\n"
  )
)

cat(
  "\n============================================================\n"
)

cat(
  "FINAL CHECKPOINTS\n"
)

cat(
  "============================================================\n"
)

print(
  checkpoints
)

cat(
  "\n============================================================\n"
)

cat(
  "SESSION INFO\n"
)

cat(
  "============================================================\n"
)

print(
  sessionInfo()
)

sink()

# Separate session info
writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    LOGDIR,
    "sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n"
)

cat(
  "STRICT-SOURCE ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "Main report:\n",
  report_file,
  "\n"
)

cat(
  "Final checkpoints:\n",
  file.path(
    TABDIR,
    "20_FINAL_CHECKPOINTS.csv"
  ),
  "\n"
)
