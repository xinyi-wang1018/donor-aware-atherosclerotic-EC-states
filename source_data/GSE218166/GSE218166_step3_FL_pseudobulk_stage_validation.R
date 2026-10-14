#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260814)

BASE <- "/root/GSE218166"

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

OUT  <- file.path(BASE, "step3_FL_pseudobulk_stage_validation")
TAB  <- file.path(OUT, "tables")
PLOT <- file.path(OUT, "plots")
OBJ  <- file.path(OUT, "objects")
LOG  <- file.path(OUT, "logs")

for (d in c(OUT, TAB, PLOT, OBJ, LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(LOG, "step3_console.log"), split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("GSE218166 STEP 3\n")
cat("FL donor-stage pseudobulk Normal -> Early -> Advanced validation\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

for (fp in c(MAIN_RDS, EXPANDED_RDS)) {
  if (!file.exists(fp)) stop("Required RDS missing: ", fp)
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("Bioconductor package 'edgeR' is required.")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Bioconductor package 'limma' is required.")
}

# -------------------------------------------------------------------------
# 0. PRE-SPECIFIED EXTERNAL-VALIDATION TARGETS
#
# These signatures are fixed BEFORE stage testing and follow the project's
# already-used panels. They are not selected based on GSE218166 results.
# -------------------------------------------------------------------------
signatures <- list(
  C7_homeostatic_arterial = c(
    "CXCL12", "EFNB2", "GJA4", "SOX17", "HEY1",
    "KLF2", "KLF4", "THBD", "EPAS1"
  ),

  C3_early_activation = c(
    "ANGPT2", "CCL14", "PGF", "LOX", "CTHRC1", "THBS1"
  ),

  C2_late_vascular_remodeling = c(
    "SULF1", "ITLN1", "MGP", "ELN", "COL8A1"
  ),

  ECM_remodeling_support = c(
    "FN1", "SPARC", "COL6A2", "IGFBP5"
  )
)

core_candidates <- c("ANGPT2", "LOX", "CTHRC1", "SULF1")
target_genes <- unique(c(unlist(signatures), core_candidates))

signature_expected <- data.frame(
  program = names(signatures),
  primary_expected_contrast = c(
    "Early-Normal_and_Advanced-Normal_negative",
    "Early-Normal_positive",
    "Advanced-Normal_and_Advanced-Early_positive",
    "Advanced-Normal_and_Advanced-Early_positive"
  ),
  stringsAsFactors = FALSE
)

candidate_expected <- data.frame(
  gene = core_candidates,
  expected_role = c(
    "early_activation",
    "early_activation",
    "early_activation",
    "late_remodeling"
  ),
  primary_expected_contrast = c(
    "Early-Normal_positive",
    "Early-Normal_positive",
    "Early-Normal_positive",
    "Advanced-Normal_and_Advanced-Early_positive"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  signature_expected,
  file.path(TAB, "00_prespecified_signature_expectations.csv"),
  row.names = FALSE
)

write.csv(
  candidate_expected,
  file.path(TAB, "00_prespecified_core_candidate_expectations.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 1. Robust raw-count extraction under Seurat v4/v5
# -------------------------------------------------------------------------
extract_counts <- function(sobj) {
  DefaultAssay(sobj) <- "RNA"
  assay <- sobj[["RNA"]]

  # First try standard single counts layer.
  ans <- tryCatch(
    GetAssayData(sobj, assay = "RNA", layer = "counts"),
    error = function(e) NULL
  )

  if (!is.null(ans) && nrow(ans) > 0 && ncol(ans) == ncol(sobj)) {
    return(ans[, colnames(sobj), drop = FALSE])
  }

  # Fall back to v4 slot API.
  ans <- tryCatch(
    GetAssayData(sobj, assay = "RNA", slot = "counts"),
    error = function(e) NULL
  )

  if (!is.null(ans) && nrow(ans) > 0 && ncol(ans) == ncol(sobj)) {
    return(ans[, colnames(sobj), drop = FALSE])
  }

  # Handle a split Seurat v5 assay containing counts.* layers.
  layer_names <- tryCatch(
    Layers(assay),
    error = function(e) character(0)
  )

  count_layers <- grep("^counts($|\\.)", layer_names, value = TRUE)

  if (length(count_layers) > 0) {
    mats <- lapply(count_layers, function(ll) {
      LayerData(assay, layer = ll)
    })
    ans <- do.call(cbind, mats)

    if (!all(colnames(sobj) %in% colnames(ans))) {
      stop("Split counts layers do not contain every object cell.")
    }

    ans <- ans[, colnames(sobj), drop = FALSE]
    return(ans)
  }

  stop("Could not retrieve RNA raw counts.")
}

# -------------------------------------------------------------------------
# 2. Donor-stage pseudobulk aggregation
# -------------------------------------------------------------------------
make_pseudobulk <- function(sobj, analysis_set, tech_keep = "FL") {

  sobj <- subset(sobj, subset = tech == tech_keep)

  if (ncol(sobj) == 0) {
    stop("No ", tech_keep, " cells in ", analysis_set)
  }

  required <- c("donor_id", "disease_stage", "tech")
  miss <- setdiff(required, colnames(sobj@meta.data))
  if (length(miss) > 0) {
    stop("Missing metadata in ", analysis_set, ": ", paste(miss, collapse = ", "))
  }

  sobj$disease_stage <- factor(
    as.character(sobj$disease_stage),
    levels = c("Normal", "Early", "Advanced"),
    ordered = TRUE
  )

  sobj$donor_stage_unit <- paste(
    sobj$donor_id,
    as.character(sobj$disease_stage),
    sep = "__"
  )

  counts <- extract_counts(sobj)

  groups <- factor(
    sobj$donor_stage_unit,
    levels = unique(sobj$donor_stage_unit)
  )

  membership <- Matrix::sparseMatrix(
    i = seq_along(groups),
    j = as.integer(groups),
    x = 1,
    dims = c(length(groups), nlevels(groups)),
    dimnames = list(colnames(sobj), levels(groups))
  )

  counts <- counts[, rownames(membership), drop = FALSE]
  pb_counts <- counts %*% membership
  pb_counts <- as.matrix(pb_counts)

  md <- sobj@meta.data %>%
    mutate(donor_stage_unit = as.character(donor_stage_unit)) %>%
    select(
      donor_stage_unit, donor_id, disease_stage, tech
    ) %>%
    distinct()

  cell_counts <- sobj@meta.data %>%
    mutate(donor_stage_unit = as.character(donor_stage_unit)) %>%
    count(
      donor_stage_unit,
      donor_id,
      disease_stage,
      tech,
      name = "n_EC_cells"
    )

  md <- md %>%
    left_join(cell_counts,
              by = c("donor_stage_unit", "donor_id",
                     "disease_stage", "tech")) %>%
    mutate(
      analysis_set = analysis_set,
      library_size = colSums(pb_counts)[donor_stage_unit]
    )

  md <- md[match(colnames(pb_counts), md$donor_stage_unit), , drop = FALSE]

  stopifnot(identical(colnames(pb_counts), md$donor_stage_unit))

  list(
    counts = pb_counts,
    meta = md,
    object = sobj
  )
}

main_all <- readRDS(MAIN_RDS)
expanded_all <- readRDS(EXPANDED_RDS)

main_pb <- make_pseudobulk(
  main_all,
  analysis_set = "MAIN_A_STRICT_EC",
  tech_keep = "FL"
)

expanded_pb <- make_pseudobulk(
  expanded_all,
  analysis_set = "EXPANDED_AB_EC",
  tech_keep = "FL"
)

write.csv(
  main_pb$counts,
  gzfile(file.path(TAB, "01_MAIN_FL_donor_stage_raw_pseudobulk_counts.csv.gz"))
)

write.csv(
  main_pb$meta,
  file.path(TAB, "02_MAIN_FL_donor_stage_metadata.csv"),
  row.names = FALSE
)

write.csv(
  expanded_pb$counts,
  gzfile(file.path(TAB, "03_EXPANDED_FL_donor_stage_raw_pseudobulk_counts.csv.gz"))
)

write.csv(
  expanded_pb$meta,
  file.path(TAB, "04_EXPANDED_FL_donor_stage_metadata.csv"),
  row.names = FALSE
)

cat("MAIN FL pseudobulk samples:", ncol(main_pb$counts), "\n")
print(table(main_pb$meta$disease_stage))
cat("\n")

# Expected 4 Normal + 5 Early + 7 Advanced = 16 units.
if (ncol(main_pb$counts) < 12) {
  stop("Unexpectedly few MAIN FL donor-stage pseudobulk samples.")
}

# -------------------------------------------------------------------------
# 3. Pseudobulk normalization and helper functions
# -------------------------------------------------------------------------
normalize_pb <- function(pb_counts) {
  y <- edgeR::DGEList(counts = pb_counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  list(y = y, logcpm = logcpm)
}

main_norm <- normalize_pb(main_pb$counts)
expanded_norm <- normalize_pb(expanded_pb$counts)

write.csv(
  main_norm$logcpm,
  gzfile(file.path(TAB, "05_MAIN_FL_TMM_logCPM.csv.gz"))
)

# Gene eligibility:
# Genome-wide: standard edgeR filterByExpr.
# Pre-specified validation genes are retained if they have at least
# 10 total pseudobulk counts and nonzero counts in >=3 donor-stage units.
build_keep <- function(y, meta, target_genes) {
  group <- factor(meta$disease_stage)
  keep_standard <- edgeR::filterByExpr(y, group = group)

  target_present <- rownames(y) %in% target_genes
  target_minimal <- (
    rowSums(y$counts) >= 10 &
      rowSums(y$counts > 0) >= 3 &
      target_present
  )

  keep_standard | target_minimal
}

# -------------------------------------------------------------------------
# 4. Differential-expression engines
#
# Normal vs Early: independent donor sets.
# Normal vs Advanced: independent donor sets.
# Early vs Advanced: PRIMARY paired analysis on the four donors
# present in both stages.
# Ordered stage trend: secondary all-sample analysis with donor blocking.
# -------------------------------------------------------------------------
run_unpaired <- function(pb_counts, meta, stage_num, stage_den,
                         target_genes, contrast_label) {

  idx <- as.character(meta$disease_stage) %in% c(stage_den, stage_num)
  counts <- pb_counts[, idx, drop = FALSE]
  md <- meta[idx, , drop = FALSE]

  md$disease_stage <- factor(
    as.character(md$disease_stage),
    levels = c(stage_den, stage_num)
  )

  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y)

  keep <- build_keep(y, md, target_genes)
  y2 <- y[keep, , keep.lib.sizes = FALSE]
  y2 <- edgeR::calcNormFactors(y2)

  design <- model.matrix(~0 + disease_stage, data = md)
  colnames(design) <- make.names(levels(md$disease_stage))

  v <- limma::voom(y2, design, plot = FALSE)
  fit <- limma::lmFit(v, design)

  contrast_expr <- paste0(
    make.names(stage_num),
    "-",
    make.names(stage_den)
  )

  cm <- limma::makeContrasts(
    contrasts = contrast_expr,
    levels = design
  )

  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2, robust = TRUE)

  tab <- limma::topTable(
    fit2,
    number = Inf,
    sort.by = "P"
  )

  tab$gene <- rownames(tab)
  tab$contrast <- contrast_label
  rownames(tab) <- NULL

  tab
}

run_paired_EA <- function(pb_counts, meta, target_genes) {

  early_donors <- meta$donor_id[
    as.character(meta$disease_stage) == "Early"
  ]
  advanced_donors <- meta$donor_id[
    as.character(meta$disease_stage) == "Advanced"
  ]

  paired_donors <- sort(intersect(early_donors, advanced_donors))

  if (length(paired_donors) < 3) {
    stop("Fewer than 3 paired Early-Advanced donors.")
  }

  idx <- meta$donor_id %in% paired_donors &
    as.character(meta$disease_stage) %in% c("Early", "Advanced")

  counts <- pb_counts[, idx, drop = FALSE]
  md <- meta[idx, , drop = FALSE]

  md$disease_stage <- factor(
    as.character(md$disease_stage),
    levels = c("Early", "Advanced")
  )
  md$donor_id <- factor(md$donor_id)

  # Ensure one Early and one Advanced unit per paired donor.
  pair_check <- table(md$donor_id, md$disease_stage)
  if (!all(pair_check == 1)) {
    stop("Paired Early-Advanced design is not exactly 1 unit/stage/donor.")
  }

  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y)

  keep <- build_keep(y, md, target_genes)
  y2 <- y[keep, , keep.lib.sizes = FALSE]
  y2 <- edgeR::calcNormFactors(y2)

  design <- model.matrix(~ donor_id + disease_stage, data = md)

  v <- limma::voom(y2, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit, robust = TRUE)

  coef_name <- grep(
    "^disease_stageAdvanced$",
    colnames(design),
    value = TRUE
  )

  if (length(coef_name) != 1) {
    stop("Could not identify paired Advanced-vs-Early coefficient.")
  }

  tab <- limma::topTable(
    fit,
    coef = coef_name,
    number = Inf,
    sort.by = "P"
  )

  tab$gene <- rownames(tab)
  tab$contrast <- "Advanced-Early_paired"
  tab$n_paired_donors <- length(paired_donors)
  rownames(tab) <- NULL

  list(
    table = tab,
    paired_donors = paired_donors,
    metadata = md
  )
}

run_stage_trend <- function(pb_counts, meta, target_genes) {

  md <- meta
  md$stage_num <- c(
    Normal = 0,
    Early = 1,
    Advanced = 2
  )[as.character(md$disease_stage)]

  y <- edgeR::DGEList(counts = pb_counts)
  y <- edgeR::calcNormFactors(y)

  keep <- build_keep(y, md, target_genes)
  y2 <- y[keep, , keep.lib.sizes = FALSE]
  y2 <- edgeR::calcNormFactors(y2)

  design <- model.matrix(~ stage_num, data = md)

  v0 <- limma::voom(y2, design, plot = FALSE)

  corfit <- tryCatch(
    limma::duplicateCorrelation(
      v0,
      design,
      block = factor(md$donor_id)
    ),
    error = function(e) NULL
  )

  use_block <- FALSE
  consensus_cor <- NA_real_

  if (!is.null(corfit) &&
      is.finite(corfit$consensus.correlation)) {
    consensus_cor <- corfit$consensus.correlation
    use_block <- TRUE
  }

  if (use_block) {
    v <- tryCatch(
      limma::voom(
        y2,
        design,
        plot = FALSE,
        block = factor(md$donor_id),
        correlation = consensus_cor
      ),
      error = function(e) v0
    )

    fit <- limma::lmFit(
      v,
      design,
      block = factor(md$donor_id),
      correlation = consensus_cor
    )
  } else {
    v <- v0
    fit <- limma::lmFit(v, design)
  }

  fit <- limma::eBayes(fit, robust = TRUE)

  tab <- limma::topTable(
    fit,
    coef = "stage_num",
    number = Inf,
    sort.by = "P"
  )

  tab$gene <- rownames(tab)
  tab$contrast <- "Ordered_stage_trend_0_1_2"
  tab$duplicateCorrelation_used <- use_block
  tab$consensus_correlation <- consensus_cor
  rownames(tab) <- NULL

  list(
    table = tab,
    consensus_correlation = consensus_cor,
    duplicateCorrelation_used = use_block
  )
}

# MAIN inference
main_EN <- run_unpaired(
  main_pb$counts, main_pb$meta,
  stage_num = "Early",
  stage_den = "Normal",
  target_genes = target_genes,
  contrast_label = "Early-Normal"
)

main_AN <- run_unpaired(
  main_pb$counts, main_pb$meta,
  stage_num = "Advanced",
  stage_den = "Normal",
  target_genes = target_genes,
  contrast_label = "Advanced-Normal"
)

main_EA <- run_paired_EA(
  main_pb$counts,
  main_pb$meta,
  target_genes
)

main_trend <- run_stage_trend(
  main_pb$counts,
  main_pb$meta,
  target_genes
)

main_all_de <- bind_rows(
  main_EN,
  main_AN,
  main_EA$table,
  main_trend$table
)

write.csv(
  main_all_de,
  gzfile(file.path(TAB, "06_MAIN_FL_all_gene_DE_results.csv.gz")),
  row.names = FALSE
)

writeLines(
  main_EA$paired_donors,
  file.path(TAB, "07_MAIN_FL_paired_Early_Advanced_donors.txt")
)

# -------------------------------------------------------------------------
# 5. Pre-specified target-gene results
# -------------------------------------------------------------------------
extract_target_table <- function(de_table, contrast_name,
                                 all_logcpm, target_genes) {

  all_genes <- intersect(target_genes, rownames(all_logcpm))

  descriptive <- data.frame(
    gene = all_genes,
    stringsAsFactors = FALSE
  )

  if (contrast_name == "Early-Normal") {
    md <- main_pb$meta
    num_stage <- "Early"
    den_stage <- "Normal"
  } else if (contrast_name == "Advanced-Normal") {
    md <- main_pb$meta
    num_stage <- "Advanced"
    den_stage <- "Normal"
  } else {
    md <- main_pb$meta
    num_stage <- NA_character_
    den_stage <- NA_character_
  }

  if (!is.na(num_stage)) {
    num_samples <- md$donor_stage_unit[
      as.character(md$disease_stage) == num_stage
    ]
    den_samples <- md$donor_stage_unit[
      as.character(md$disease_stage) == den_stage
    ]

    descriptive$mean_logCPM_numerator <- rowMeans(
      all_logcpm[all_genes, num_samples, drop = FALSE]
    )
    descriptive$mean_logCPM_denominator <- rowMeans(
      all_logcpm[all_genes, den_samples, drop = FALSE]
    )
    descriptive$descriptive_delta_logCPM <-
      descriptive$mean_logCPM_numerator -
      descriptive$mean_logCPM_denominator
  }

  stat <- de_table %>%
    filter(gene %in% all_genes) %>%
    select(
      gene,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B
    )

  out <- descriptive %>%
    left_join(stat, by = "gene") %>%
    mutate(
      contrast = contrast_name,
      tested_by_limma = !is.na(P.Value)
    )

  out
}

target_EN <- extract_target_table(
  main_EN,
  "Early-Normal",
  main_norm$logcpm,
  target_genes
)

target_AN <- extract_target_table(
  main_AN,
  "Advanced-Normal",
  main_norm$logcpm,
  target_genes
)

# Paired EA descriptive effects for ALL targets.
paired_donors <- main_EA$paired_donors
paired_md <- main_EA$metadata

ea_descriptive <- lapply(
  intersect(target_genes, rownames(main_norm$logcpm)),
  function(g) {
    diffs <- sapply(paired_donors, function(dd) {
      e_sample <- paired_md$donor_stage_unit[
        paired_md$donor_id == dd &
          as.character(paired_md$disease_stage) == "Early"
      ]
      a_sample <- paired_md$donor_stage_unit[
        paired_md$donor_id == dd &
          as.character(paired_md$disease_stage) == "Advanced"
      ]

      main_norm$logcpm[g, a_sample] -
        main_norm$logcpm[g, e_sample]
    })

    data.frame(
      gene = g,
      mean_paired_delta_logCPM = mean(diffs),
      median_paired_delta_logCPM = median(diffs),
      n_positive_pairs = sum(diffs > 0),
      n_negative_pairs = sum(diffs < 0),
      n_pairs = length(diffs),
      stringsAsFactors = FALSE
    )
  }
) %>% bind_rows()

target_EA <- main_EA$table %>%
  filter(gene %in% target_genes) %>%
  select(
    gene,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B,
    n_paired_donors
  ) %>%
  full_join(ea_descriptive, by = "gene") %>%
  mutate(
    contrast = "Advanced-Early_paired",
    tested_by_limma = !is.na(P.Value)
  )

target_trend <- main_trend$table %>%
  filter(gene %in% target_genes) %>%
  select(
    gene,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B,
    duplicateCorrelation_used,
    consensus_correlation
  ) %>%
  mutate(
    contrast = "Ordered_stage_trend_0_1_2",
    tested_by_limma = TRUE
  )

main_target_all <- bind_rows(
  target_EN,
  target_AN,
  target_EA,
  target_trend
)

# Target-panel FDR: correction ONLY across the pre-specified genes
# within each contrast, separately from genome-wide FDR.
main_target_all <- main_target_all %>%
  group_by(contrast) %>%
  mutate(
    target_panel_FDR = p.adjust(P.Value, method = "BH")
  ) %>%
  ungroup()

write.csv(
  main_target_all,
  file.path(TAB, "08_MAIN_FL_prespecified_target_gene_results.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 6. Program scores at DONOR-STAGE level
#
# Gene-wise z-standardization is done across pseudobulk samples, then
# program score = mean z across pre-specified program genes.
# -------------------------------------------------------------------------
make_program_scores <- function(logcpm, meta, signatures, analysis_set) {

  gene_z <- t(scale(t(logcpm)))
  gene_z[!is.finite(gene_z)] <- NA_real_

  out <- lapply(names(signatures), function(nm) {
    genes <- intersect(signatures[[nm]], rownames(gene_z))

    if (length(genes) == 0) return(NULL)

    vals <- colMeans(
      gene_z[genes, , drop = FALSE],
      na.rm = TRUE
    )

    data.frame(
      donor_stage_unit = colnames(logcpm),
      program = nm,
      score = as.numeric(vals),
      n_genes_available = length(genes),
      genes_available = paste(genes, collapse = ";"),
      analysis_set = analysis_set,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  out %>%
    left_join(
      meta %>%
        select(
          donor_stage_unit,
          donor_id,
          disease_stage,
          tech,
          n_EC_cells
        ),
      by = "donor_stage_unit"
    )
}

main_program_scores <- make_program_scores(
  main_norm$logcpm,
  main_pb$meta,
  signatures,
  "MAIN_A_STRICT_EC"
)

expanded_program_scores <- make_program_scores(
  expanded_norm$logcpm,
  expanded_pb$meta,
  signatures,
  "EXPANDED_AB_EC"
)

write.csv(
  main_program_scores,
  file.path(TAB, "09_MAIN_FL_donor_stage_program_scores.csv"),
  row.names = FALSE
)

write.csv(
  expanded_program_scores,
  file.path(TAB, "10_EXPANDED_FL_donor_stage_program_scores.csv"),
  row.names = FALSE
)

program_pairwise_stats <- function(score_df, analysis_set) {

  programs <- unique(score_df$program)

  all_stats <- lapply(programs, function(pg) {
    d <- score_df %>% filter(program == pg)

    get_vals <- function(stage) {
      d$score[as.character(d$disease_stage) == stage]
    }

    N <- get_vals("Normal")
    E <- get_vals("Early")
    A <- get_vals("Advanced")

    # Welch tests are secondary; direction/effect size is reported regardless.
    en_t <- tryCatch(t.test(E, N), error = function(e) NULL)
    an_t <- tryCatch(t.test(A, N), error = function(e) NULL)

    # Paired Early-Advanced on repeated donors.
    e_df <- d %>%
      filter(as.character(disease_stage) == "Early") %>%
      select(donor_id, E = score)

    a_df <- d %>%
      filter(as.character(disease_stage) == "Advanced") %>%
      select(donor_id, A = score)

    pair <- inner_join(e_df, a_df, by = "donor_id")

    ea_t <- if (nrow(pair) >= 3) {
      tryCatch(
        t.test(pair$A, pair$E, paired = TRUE),
        error = function(e) NULL
      )
    } else NULL

    bind_rows(
      data.frame(
        analysis_set = analysis_set,
        program = pg,
        contrast = "Early-Normal",
        effect = mean(E) - mean(N),
        p_value = ifelse(is.null(en_t), NA, en_t$p.value),
        n_num = length(E),
        n_den = length(N),
        n_paired = NA_integer_,
        n_pairs_expected_positive = NA_integer_
      ),
      data.frame(
        analysis_set = analysis_set,
        program = pg,
        contrast = "Advanced-Normal",
        effect = mean(A) - mean(N),
        p_value = ifelse(is.null(an_t), NA, an_t$p.value),
        n_num = length(A),
        n_den = length(N),
        n_paired = NA_integer_,
        n_pairs_expected_positive = NA_integer_
      ),
      data.frame(
        analysis_set = analysis_set,
        program = pg,
        contrast = "Advanced-Early_paired",
        effect = ifelse(
          nrow(pair) > 0,
          mean(pair$A - pair$E),
          NA_real_
        ),
        p_value = ifelse(is.null(ea_t), NA, ea_t$p.value),
        n_num = NA_integer_,
        n_den = NA_integer_,
        n_paired = nrow(pair),
        n_pairs_expected_positive = ifelse(
          nrow(pair) > 0,
          sum(pair$A - pair$E > 0),
          NA_integer_
        )
      )
    )
  }) %>% bind_rows()

  all_stats %>%
    group_by(contrast) %>%
    mutate(program_FDR = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

main_program_stats <- program_pairwise_stats(
  main_program_scores,
  "MAIN_A_STRICT_EC"
)

expanded_program_stats <- program_pairwise_stats(
  expanded_program_scores,
  "EXPANDED_AB_EC"
)

program_stats_all <- bind_rows(
  main_program_stats,
  expanded_program_stats
)

write.csv(
  program_stats_all,
  file.path(TAB, "11_program_stage_comparison_results.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 7. Expanded A+B sensitivity for TARGET GENES
# -------------------------------------------------------------------------
expanded_EN <- run_unpaired(
  expanded_pb$counts, expanded_pb$meta,
  stage_num = "Early",
  stage_den = "Normal",
  target_genes = target_genes,
  contrast_label = "Early-Normal"
)

expanded_AN <- run_unpaired(
  expanded_pb$counts, expanded_pb$meta,
  stage_num = "Advanced",
  stage_den = "Normal",
  target_genes = target_genes,
  contrast_label = "Advanced-Normal"
)

expanded_EA <- run_paired_EA(
  expanded_pb$counts,
  expanded_pb$meta,
  target_genes
)

expanded_target <- bind_rows(
  expanded_EN %>%
    filter(gene %in% target_genes),
  expanded_AN %>%
    filter(gene %in% target_genes),
  expanded_EA$table %>%
    filter(gene %in% target_genes)
) %>%
  group_by(contrast) %>%
  mutate(target_panel_FDR = p.adjust(P.Value, method = "BH")) %>%
  ungroup() %>%
  mutate(analysis_set = "EXPANDED_AB_EC")

write.csv(
  expanded_target,
  file.path(TAB, "12_EXPANDED_FL_prespecified_target_gene_results.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8. Main-vs-expanded concordance for core candidates
# -------------------------------------------------------------------------
main_core <- main_target_all %>%
  filter(
    gene %in% core_candidates,
    contrast %in%
      c("Early-Normal", "Advanced-Normal", "Advanced-Early_paired")
  ) %>%
  select(
    gene,
    contrast,
    main_logFC = logFC,
    main_P = P.Value,
    main_target_FDR = target_panel_FDR
  )

expanded_core <- expanded_target %>%
  filter(gene %in% core_candidates) %>%
  select(
    gene,
    contrast,
    expanded_logFC = logFC,
    expanded_P = P.Value,
    expanded_target_FDR = target_panel_FDR
  )

core_concordance <- full_join(
  main_core,
  expanded_core,
  by = c("gene", "contrast")
) %>%
  mutate(
    same_direction =
      sign(main_logFC) == sign(expanded_logFC)
  )

write.csv(
  core_concordance,
  file.path(TAB, "13_core_candidate_MAIN_vs_EXPANDED_concordance.csv"),
  row.names = FALSE
)

program_concordance <- main_program_stats %>%
  select(
    program,
    contrast,
    main_effect = effect,
    main_P = p_value,
    main_FDR = program_FDR
  ) %>%
  left_join(
    expanded_program_stats %>%
      select(
        program,
        contrast,
        expanded_effect = effect,
        expanded_P = p_value,
        expanded_FDR = program_FDR
      ),
    by = c("program", "contrast")
  ) %>%
  mutate(
    same_direction =
      sign(main_effect) == sign(expanded_effect)
  )

write.csv(
  program_concordance,
  file.path(TAB, "14_program_MAIN_vs_EXPANDED_concordance.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 9. Evidence classification WITHOUT forcing significance
# -------------------------------------------------------------------------
get_prog_effect <- function(program, contrast, set = "MAIN_A_STRICT_EC") {
  x <- program_stats_all %>%
    filter(
      analysis_set == set,
      .data$program == program,
      .data$contrast == contrast
    )
  if (nrow(x) == 0) return(NA_real_)
  x$effect[1]
}

get_prog_p <- function(program, contrast, set = "MAIN_A_STRICT_EC") {
  x <- program_stats_all %>%
    filter(
      analysis_set == set,
      .data$program == program,
      .data$contrast == contrast
    )
  if (nrow(x) == 0) return(NA_real_)
  x$p_value[1]
}

program_evidence <- data.frame(
  program = names(signatures),
  stringsAsFactors = FALSE
)

program_evidence$main_Early_minus_Normal <- sapply(
  program_evidence$program,
  get_prog_effect,
  contrast = "Early-Normal"
)

program_evidence$main_Advanced_minus_Normal <- sapply(
  program_evidence$program,
  get_prog_effect,
  contrast = "Advanced-Normal"
)

program_evidence$main_Advanced_minus_Early_paired <- sapply(
  program_evidence$program,
  get_prog_effect,
  contrast = "Advanced-Early_paired"
)

program_evidence$expanded_Early_minus_Normal <- sapply(
  program_evidence$program,
  get_prog_effect,
  contrast = "Early-Normal",
  set = "EXPANDED_AB_EC"
)

program_evidence$expanded_Advanced_minus_Normal <- sapply(
  program_evidence$program,
  get_prog_effect,
  contrast = "Advanced-Normal",
  set = "EXPANDED_AB_EC"
)

program_evidence$expanded_Advanced_minus_Early_paired <- sapply(
  program_evidence$program,
  get_prog_effect,
  contrast = "Advanced-Early_paired",
  set = "EXPANDED_AB_EC"
)

program_evidence$directional_evidence <- sapply(
  program_evidence$program,
  function(pg) {
    r <- program_evidence[
      program_evidence$program == pg, ,
      drop = FALSE
    ]

    if (pg == "C7_homeostatic_arterial") {
      main_ok <- (
        r$main_Early_minus_Normal < 0 &&
          r$main_Advanced_minus_Normal < 0
      )
      sens_ok <- (
        r$expanded_Early_minus_Normal < 0 &&
          r$expanded_Advanced_minus_Normal < 0
      )
    } else if (pg == "C3_early_activation") {
      main_ok <- r$main_Early_minus_Normal > 0
      sens_ok <- r$expanded_Early_minus_Normal > 0
    } else {
      main_ok <- (
        r$main_Advanced_minus_Normal > 0 &&
          r$main_Advanced_minus_Early_paired > 0
      )
      sens_ok <- (
        r$expanded_Advanced_minus_Normal > 0 &&
          r$expanded_Advanced_minus_Early_paired > 0
      )
    }

    if (isTRUE(main_ok) && isTRUE(sens_ok)) {
      "directionally_reproduced_main_and_expanded"
    } else if (isTRUE(main_ok)) {
      "directionally_supported_in_main_only"
    } else {
      "expected_direction_not_reproduced_in_main"
    }
  }
)

write.csv(
  program_evidence,
  file.path(TAB, "15_primary_program_evidence_summary.csv"),
  row.names = FALSE
)

# Candidate-level evidence
candidate_summary <- core_concordance %>%
  filter(
    (gene %in% c("ANGPT2","LOX","CTHRC1") &
       contrast == "Early-Normal") |
      (gene == "SULF1" &
         contrast %in%
           c("Advanced-Normal","Advanced-Early_paired"))
  ) %>%
  left_join(candidate_expected, by = "gene") %>%
  mutate(
    expected_direction_met_main = case_when(
      gene %in% c("ANGPT2","LOX","CTHRC1") ~ main_logFC > 0,
      gene == "SULF1" ~ main_logFC > 0,
      TRUE ~ NA
    ),
    expected_direction_met_expanded = case_when(
      gene %in% c("ANGPT2","LOX","CTHRC1") ~ expanded_logFC > 0,
      gene == "SULF1" ~ expanded_logFC > 0,
      TRUE ~ NA
    ),
    evidence_label = case_when(
      expected_direction_met_main &
        expected_direction_met_expanded &
        !is.na(main_P) &
        main_P < 0.05
        ~ "expected_direction_main_and_sensitivity_with_nominal_main_support",

      expected_direction_met_main &
        expected_direction_met_expanded
        ~ "expected_direction_main_and_sensitivity_without_nominal_main_support",

      expected_direction_met_main
        ~ "expected_direction_main_only",

      TRUE
        ~ "expected_direction_not_supported"
    )
  )

write.csv(
  candidate_summary,
  file.path(TAB, "16_core_candidate_evidence_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 10. 3DE: DESCRIPTIVE ONLY because Early has one donor after final freeze.
# -------------------------------------------------------------------------
make_3de_descriptive <- function(sobj, signatures) {
  pb <- make_pseudobulk(
    sobj,
    analysis_set = "MAIN_A_STRICT_EC_3DE",
    tech_keep = "3DE"
  )
  norm <- normalize_pb(pb$counts)
  scores <- make_program_scores(
    norm$logcpm,
    pb$meta,
    signatures,
    "MAIN_A_STRICT_EC_3DE"
  )

  stage_summary <- scores %>%
    group_by(program, disease_stage) %>%
    summarise(
      n_donors = n_distinct(donor_id),
      mean_score = mean(score),
      median_score = median(score),
      .groups = "drop"
    )

  list(pb = pb, scores = scores, stage_summary = stage_summary)
}

de3 <- make_3de_descriptive(main_all, signatures)

write.csv(
  de3$pb$meta,
  file.path(TAB, "17_3DE_MAIN_donor_stage_metadata_DESCRIPTIVE_ONLY.csv"),
  row.names = FALSE
)

write.csv(
  de3$stage_summary,
  file.path(TAB, "18_3DE_program_stage_DESCRIPTIVE_ONLY.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 11. Figures - ALL PURE WHITE BACKGROUND
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

# 11a. Pseudobulk PCA
pca_input <- main_norm$logcpm
rv <- apply(pca_input, 1, var, na.rm = TRUE)
top_var <- names(sort(rv, decreasing = TRUE))[
  seq_len(min(1000, sum(is.finite(rv) & rv > 0)))
]

pca <- prcomp(
  t(pca_input[top_var, , drop = FALSE]),
  center = TRUE,
  scale. = FALSE
)

pca_df <- data.frame(
  donor_stage_unit = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  stringsAsFactors = FALSE
) %>%
  left_join(
    main_pb$meta %>%
      select(donor_stage_unit, donor_id, disease_stage),
    by = "donor_stage_unit"
  )

p1 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    shape = disease_stage,
    label = donor_id
  )
) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  theme_classic(base_size = 11) +
  labs(
    title = "FL main strict EC donor-stage pseudobulk PCA",
    subtitle = "Biological replicate = donor × disease stage"
  )

save_white(p1, "01_MAIN_FL_pseudobulk_PCA.png", 10, 7)

# 11b. Program scores
p2 <- ggplot(
  main_program_scores,
  aes(
    x = disease_stage,
    y = score,
    group = disease_stage
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2) +
  facet_wrap(~program, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 10) +
  labs(
    title = "Pre-specified EC program scores across disease stages",
    x = NULL,
    y = "Pseudobulk program z-score"
  )

save_white(p2, "02_MAIN_FL_program_scores_by_stage.png", 12, 8)

# 11c. Core candidate expression
core_long <- as.data.frame(
  t(
    main_norm$logcpm[
      intersect(core_candidates, rownames(main_norm$logcpm)),
      ,
      drop = FALSE
    ]
  )
)

core_long$donor_stage_unit <- rownames(core_long)

core_long <- core_long %>%
  pivot_longer(
    cols = all_of(intersect(core_candidates, colnames(core_long))),
    names_to = "gene",
    values_to = "logCPM"
  ) %>%
  left_join(
    main_pb$meta %>%
      select(donor_stage_unit, donor_id, disease_stage),
    by = "donor_stage_unit"
  )

p3 <- ggplot(
  core_long,
  aes(
    x = disease_stage,
    y = logCPM
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2) +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 10) +
  labs(
    title = "Core candidate genes in FL strict-EC pseudobulk",
    x = NULL,
    y = "TMM logCPM"
  )

save_white(p3, "03_MAIN_FL_core_candidates_by_stage.png", 11, 8)

# 11d. Target-gene heatmap
target_heat_genes <- intersect(
  target_genes,
  rownames(main_norm$logcpm)
)

heat <- main_norm$logcpm[
  target_heat_genes,
  ,
  drop = FALSE
]

heat_z <- t(scale(t(heat)))
heat_z[!is.finite(heat_z)] <- 0

heat_df <- as.data.frame(as.table(heat_z))
colnames(heat_df) <- c("gene", "donor_stage_unit", "z")

heat_df <- heat_df %>%
  left_join(
    main_pb$meta %>%
      select(donor_stage_unit, disease_stage),
    by = "donor_stage_unit"
  )

sample_order <- main_pb$meta %>%
  arrange(disease_stage, donor_id) %>%
  pull(donor_stage_unit)

heat_df$donor_stage_unit <- factor(
  heat_df$donor_stage_unit,
  levels = sample_order
)

p4 <- ggplot(
  heat_df,
  aes(
    x = donor_stage_unit,
    y = gene,
    fill = z
  )
) +
  geom_tile() +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "Pre-specified validation genes across donor-stage pseudobulks",
    fill = "Gene z"
  )

save_white(p4, "04_MAIN_FL_target_gene_heatmap.png", 14, 10)

# 11e. Core candidate effect plot
core_effect_plot <- core_concordance %>%
  filter(
    contrast %in%
      c("Early-Normal","Advanced-Normal","Advanced-Early_paired")
  ) %>%
  pivot_longer(
    cols = c(main_logFC, expanded_logFC),
    names_to = "analysis",
    values_to = "logFC"
  )

p5 <- ggplot(
  core_effect_plot,
  aes(
    x = contrast,
    y = logFC,
    shape = analysis
  )
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 2.8
  ) +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    title = "Core candidate pseudobulk effects: strict vs expanded EC",
    x = NULL,
    y = "limma logFC"
  )

save_white(p5, "05_core_candidate_effects_MAIN_vs_EXPANDED.png", 12, 8)

# 11f. Paired Early -> Advanced trajectories for programs
paired_program <- main_program_scores %>%
  filter(
    donor_id %in% main_EA$paired_donors,
    as.character(disease_stage) %in% c("Early", "Advanced")
  )

p6 <- ggplot(
  paired_program,
  aes(
    x = disease_stage,
    y = score,
    group = donor_id
  )
) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~program, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 10) +
  labs(
    title = "Matched-donor Early to Advanced program changes",
    subtitle = paste(
      "Paired donors:",
      paste(main_EA$paired_donors, collapse = ", ")
    ),
    x = NULL,
    y = "Program z-score"
  )

save_white(p6, "06_paired_Early_Advanced_program_changes.png", 12, 8)

# -------------------------------------------------------------------------
# 12. Compact final decision report
# -------------------------------------------------------------------------
candidate_lines <- candidate_summary %>%
  mutate(
    line = paste0(
      gene, " | ", contrast,
      " | main logFC=", sprintf("%.3f", main_logFC),
      " | main P=", ifelse(
        is.na(main_P), "NA",
        format(main_P, digits = 3, scientific = TRUE)
      ),
      " | expanded logFC=", sprintf("%.3f", expanded_logFC),
      " | ", evidence_label
    )
  ) %>%
  pull(line)

program_lines <- program_evidence %>%
  mutate(
    line = paste0(
      program,
      " | E-N=", sprintf("%.3f", main_Early_minus_Normal),
      " | A-N=", sprintf("%.3f", main_Advanced_minus_Normal),
      " | A-E paired=", sprintf("%.3f", main_Advanced_minus_Early_paired),
      " | ", directional_evidence
    )
  ) %>%
  pull(line)

summary_lines <- c(
  "============================================================",
  "GSE218166 STEP3 FL PSEUDOBULK STAGE VALIDATION",
  "============================================================",
  paste0("Run time: ", format(Sys.time())),
  "",
  "PRIMARY BIOLOGICAL REPLICATE",
  "donor x disease-stage pseudobulk; individual cells are not treated as independent replicates.",
  "",
  paste0("MAIN FL pseudobulk units: ", ncol(main_pb$counts)),
  paste0(
    "Normal / Early / Advanced units: ",
    paste(
      as.integer(table(main_pb$meta$disease_stage)),
      collapse = " / "
    )
  ),
  paste0(
    "Paired Early-Advanced donors: ",
    paste(main_EA$paired_donors, collapse = ", ")
  ),
  "",
  "PRIMARY CONTRASTS",
  "Early-Normal: unpaired donor-level pseudobulk.",
  "Advanced-Normal: unpaired donor-level pseudobulk.",
  "Advanced-Early: paired analysis restricted to donors measured at both stages.",
  "Ordered 0-1-2 trend: secondary analysis; duplicateCorrelation donor blocking used when estimable.",
  "",
  paste0(
    "Ordered-trend donor correlation estimate: ",
    ifelse(
      is.na(main_trend$consensus_correlation),
      "NA",
      sprintf("%.3f", main_trend$consensus_correlation)
    )
  ),
  "",
  "PROGRAM DIRECTIONAL SUMMARY",
  program_lines,
  "",
  "CORE CANDIDATE SUMMARY",
  candidate_lines,
  "",
  "INTERPRETATION RULE",
  "Nominal/FDR significance is reported, but evidence is not called positive solely because P<0.05.",
  "Direction must agree with the pre-specified biological expectation, and MAIN-versus-EXPANDED sensitivity concordance is reported.",
  "3DE Early remains descriptive only because it has one donor after final EC freeze.",
  "",
  "NEXT",
  "Upload GSE218166_step3_REVIEW_SMALL.tar.gz for biological interpretation and decision on the next validation step."
)

writeLines(
  summary_lines,
  file.path(OUT, "Step3_DECISION_SUMMARY.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, "sessionInfo.txt")
)

# Save a compact R object with pseudobulk products for later Step3B.
saveRDS(
  list(
    main_counts = main_pb$counts,
    main_meta = main_pb$meta,
    main_logCPM = main_norm$logcpm,
    expanded_counts = expanded_pb$counts,
    expanded_meta = expanded_pb$meta,
    expanded_logCPM = expanded_norm$logcpm,
    signatures = signatures,
    core_candidates = core_candidates
  ),
  file.path(OBJ, "GSE218166_step3_pseudobulk_products.rds"),
  compress = "gzip"
)

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\nCompleted:", format(Sys.time()), "\n")
