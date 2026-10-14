#!/usr/bin/env bash
set -euo pipefail

OUT="/root/AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919"
R="/root/run_C3_contamination_audit_CORRECTED_20260919.R"
BUNDLE="/root/AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_SMALL_20260919.tar.gz"

rm -rf "$OUT" "$BUNDLE"
mkdir -p "$OUT"

cat > "$R" <<'RSCRIPT'
options(stringsAsFactors = FALSE)
set.seed(20260919)

OUT <- "/root/AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

zz <- file(file.path(OUT, "RUN_LOG.txt"), open="wt")
sink(zz, type="output", split=TRUE)
sink(zz, type="message")
on.exit({
  sink(type="message")
  sink(type="output")
  close(zz)
}, add=TRUE)

cat("CORRECTED C3 contamination / identity audit\n")
cat("Time:", as.character(Sys.time()), "\n\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
})

EC_RDS <- "/root/AS_EC_state_scoring/AS_EC_highconfidence_state_scored.rds"
DONOR_MAP <- "/root/AS_EC_19day_sprint/04_donor_program_validation/04_sample_to_donor_mapping.csv"

if (!file.exists(EC_RDS)) stop("Frozen EC object not found: ", EC_RDS)
if (!file.exists(DONOR_MAP)) stop("Frozen donor map not found: ", DONOR_MAP)

STRICT_DATASETS <- c("GSE155512","GSE253904")
CORE_STATES <- c("C1","C2","C3","C4","C5","C7")
EXPECTED <- c(strict_highconf=6152L, strict_core=4435L,
              C1=1440L, C2=1181L, C3=696L, C4=624L, C5=324L, C7=170L)

write_csv <- function(x, path) {
  write.csv(as.data.frame(x), path, row.names=FALSE, quote=TRUE, na="")
  cat("WRITE:", path, "[", nrow(as.data.frame(x)), "rows]\n")
}

first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% colnames(x)]
  if (length(hit)) hit[[1]] else NA_character_
}

canonical_state <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  for (st in CORE_STATES) {
    idx <- grepl(paste0("(^|[^0-9])", st, "([^0-9]|$)"), x, ignore.case=TRUE) |
      toupper(x) == st |
      sub("^[Cc]", "", x) == sub("^C", "", st)
    out[idx] <- st
  }
  out[is.na(out) & grepl("ROOT.*C7|C7.*ROOT|HOMEOST", x, ignore.case=TRUE)] <- "C7"
  out[is.na(out) & grepl("EARLY.*C3|C3.*EARLY|ACTIV", x, ignore.case=TRUE)] <- "C3"
  out[is.na(out) & grepl("TRANSITION.*C5|C5.*TRANSITION", x, ignore.case=TRUE)] <- "C5"
  out
}

attach_donor <- function(md) {
  dm <- read.csv(DONOR_MAP, check.names=FALSE, stringsAsFactors=FALSE)
  req <- c("dataset_id","sample_id","donor_id")
  if (!all(req %in% names(dm))) stop("Donor map missing: ", paste(setdiff(req,names(dm)),collapse=", "))
  dm <- unique(dm[,req,drop=FALSE])
  dkey <- paste(dm$dataset_id, dm$sample_id, sep="|||MAP|||")
  ckey <- paste(md$dataset_id, md$sample_id, sep="|||MAP|||")
  md$donor_id <- dm$donor_id[match(ckey,dkey)]
  if (anyNA(md$donor_id)) stop("Some strict EC cells could not be mapped to donor.")
  md
}

get_layer <- function(obj, layer) {
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) {
    obj <- tryCatch(SeuratObject::JoinLayers(obj, assay="RNA"), error=function(e) obj)
  }
  x <- tryCatch(SeuratObject::LayerData(obj, assay="RNA", layer=layer), error=function(e) NULL)
  if (is.null(x) || !nrow(x)) {
    x <- tryCatch(Seurat::GetAssayData(obj, assay="RNA", slot=layer), error=function(e) NULL)
  }
  if (is.null(x) || !nrow(x)) stop("RNA layer unavailable: ", layer)
  x
}

# ------------------ LOAD EXACT FROZEN EC OBJECT ------------------
obj <- readRDS(EC_RDS)
md0 <- obj@meta.data
if (!all(c("dataset_id","sample_id") %in% names(md0))) stop("EC object lacks dataset_id/sample_id.")

state_col <- first_existing(md0, c("EC_cluster","EC_cluster_std","paper_state","EC_refined_state"))
if (is.na(state_col)) stop("Frozen EC state column not found.")
cat("Frozen EC object:", EC_RDS, "\n")
cat("State source column:", state_col, "\n")

state_all <- canonical_state(md0[[state_col]])
names(state_all) <- rownames(md0)

strict_cells <- rownames(md0)[as.character(md0$dataset_id) %in% STRICT_DATASETS]
strict <- subset(obj, cells=strict_cells)
md <- strict@meta.data
md$paper_state <- state_all[rownames(md)]
md <- attach_donor(md)
strict@meta.data <- md

state_tab <- table(factor(md$paper_state[md$paper_state %in% CORE_STATES], levels=CORE_STATES))
qc <- data.frame(
  metric=c("strict_highconf","strict_core",CORE_STATES),
  observed=c(nrow(md), sum(md$paper_state %in% CORE_STATES), as.integer(state_tab)),
  expected=c(EXPECTED["strict_highconf"], EXPECTED["strict_core"], EXPECTED[CORE_STATES])
)
qc$PASS <- qc$observed == qc$expected
write_csv(qc, file.path(OUT,"00_HARD_QC_state_counts.csv"))

cat("\nHARD QC:\n")
print(qc)
if (!all(qc$PASS)) {
  stop("HARD QC FAILED. Refusing contamination interpretation.")
}

core_cells <- rownames(md)[md$paper_state %in% CORE_STATES]
core <- subset(strict, cells=core_cells)
cmd <- core@meta.data
cmd$paper_state <- md[rownames(cmd),"paper_state"]
cmd$donor_id <- md[rownames(cmd),"donor_id"]
core@meta.data <- cmd

# ------------------ EXPRESSION + PROGRAM SCORES ------------------
dat <- get_layer(core, "data")
counts <- get_layer(core, "counts")
cmd <- core@meta.data

gene_sets <- list(
  EC_identity = c("PECAM1","VWF","CDH5","KDR","ESAM","EMCN","ENG","RAMP2","ROBO4","ECSCR","ERG","EPAS1","KLF2","THBD"),
  Fibroblast = c("COL1A1","COL1A2","COL3A1","DCN","LUM","COL6A1","COL6A2","COL6A3","PDGFRA","FAP","COL5A1","COL5A2"),
  SMC = c("ACTA2","TAGLN","MYH11","CNN1","MYL9","TPM2","CALD1"),
  Pericyte = c("RGS5","CSPG4","PDGFRB","MCAM","NOTCH3","DES"),
  Myeloid = c("PTPRC","LYZ","LST1","TYROBP","FCER1G","CTSS")
)

# Gene-wise z-score across strict-core EC, then mean z score per program.
score_program <- function(gs) {
  g <- intersect(gs, rownames(dat))
  if (length(g) < 2) return(rep(NA_real_, ncol(dat)))
  mat <- as.matrix(dat[g,,drop=FALSE])
  mu <- rowMeans(mat)
  sdv <- apply(mat,1,sd)
  keep <- is.finite(sdv) & sdv > 0
  if (sum(keep) < 2) return(rep(NA_real_, ncol(mat)))
  z <- sweep(mat[keep,,drop=FALSE],1,mu[keep],"-")
  z <- sweep(z,1,sdv[keep],"/")
  colMeans(z)
}

for (nm in names(gene_sets)) {
  cmd[[paste0(nm,"_z")]] <- score_program(gene_sets[[nm]])
}
core@meta.data <- cmd

# Marker detection count/positivity for transparent identity audit.
marker_genes <- unique(unlist(gene_sets))
marker_genes <- intersect(marker_genes, rownames(dat))
marker_rows <- list()
for (st in CORE_STATES) {
  cells <- rownames(cmd)[cmd$paper_state==st]
  for (g in marker_genes) {
    x <- as.numeric(dat[g,cells])
    marker_rows[[length(marker_rows)+1]] <- data.frame(
      state=st, gene=g,
      program=names(gene_sets)[vapply(gene_sets,function(s) g %in% s,logical(1))][1],
      n_cells=length(cells),
      pct_expr_gt0=100*mean(x>0),
      mean_logexpr=mean(x),
      median_logexpr=median(x)
    )
  }
}
marker_df <- do.call(rbind, marker_rows)
write_csv(marker_df, file.path(OUT,"01_marker_expression_by_state.csv"))

score_cols <- paste0(names(gene_sets),"_z")
score_summary <- do.call(rbind,lapply(CORE_STATES,function(st){
  z <- cmd[cmd$paper_state==st,,drop=FALSE]
  out <- data.frame(state=st,n=nrow(z))
  for (cc in score_cols) {
    out[[paste0(cc,"_mean")]] <- mean(z[[cc]],na.rm=TRUE)
    out[[paste0(cc,"_median")]] <- median(z[[cc]],na.rm=TRUE)
  }
  out
}))
write_csv(score_summary, file.path(OUT,"02_program_scores_by_state.csv"))

# ------------------ FIXED C7-REFERENCED RISK THRESHOLDS ------------------
c7 <- cmd[cmd$paper_state=="C7",,drop=FALSE]
qfun <- function(x,p) as.numeric(quantile(x,p,na.rm=TRUE,type=7))

thr <- data.frame(
  metric=c("EC_identity_q25","Fibroblast_q95","SMC_q95","Pericyte_q95","Myeloid_q95"),
  threshold=c(qfun(c7$EC_identity_z,.25),
              qfun(c7$Fibroblast_z,.95),
              qfun(c7$SMC_z,.95),
              qfun(c7$Pericyte_z,.95),
              qfun(c7$Myeloid_z,.95))
)
write_csv(thr, file.path(OUT,"03_C7_reference_thresholds.csv"))
TV <- setNames(thr$threshold,thr$metric)

cmd$Fibroblast_high <- cmd$Fibroblast_z > TV["Fibroblast_q95"]
cmd$SMC_high <- cmd$SMC_z > TV["SMC_q95"]
cmd$Pericyte_high <- cmd$Pericyte_z > TV["Pericyte_q95"]
cmd$Myeloid_high <- cmd$Myeloid_z > TV["Myeloid_q95"]
cmd$any_nonEC_high <- cmd$Fibroblast_high | cmd$SMC_high | cmd$Pericyte_high | cmd$Myeloid_high
cmd$EC_low <- cmd$EC_identity_z < TV["EC_identity_q25"]
# Main suspicious flag requires both non-EC-high and relatively low EC identity.
cmd$combined_nonEC_risk <- cmd$any_nonEC_high & cmd$EC_low
core@meta.data <- cmd

risk_summary <- do.call(rbind,lapply(CORE_STATES,function(st){
  z <- cmd[cmd$paper_state==st,,drop=FALSE]
  data.frame(
    state=st,n=nrow(z),
    fibroblast_high_n=sum(z$Fibroblast_high), fibroblast_high_pct=100*mean(z$Fibroblast_high),
    smc_high_n=sum(z$SMC_high), smc_high_pct=100*mean(z$SMC_high),
    pericyte_high_n=sum(z$Pericyte_high), pericyte_high_pct=100*mean(z$Pericyte_high),
    myeloid_high_n=sum(z$Myeloid_high), myeloid_high_pct=100*mean(z$Myeloid_high),
    any_nonEC_high_n=sum(z$any_nonEC_high), any_nonEC_high_pct=100*mean(z$any_nonEC_high),
    EC_low_n=sum(z$EC_low), EC_low_pct=100*mean(z$EC_low),
    combined_nonEC_risk_n=sum(z$combined_nonEC_risk),
    combined_nonEC_risk_pct=100*mean(z$combined_nonEC_risk)
  )
}))
write_csv(risk_summary, file.path(OUT,"04_nonEC_risk_by_state.csv"))

# Donor-level C3 audit.
c3md <- cmd[cmd$paper_state=="C3",,drop=FALSE]
donor_summary <- do.call(rbind,lapply(split(c3md,c3md$donor_id),function(z){
  data.frame(
    donor_id=z$donor_id[1],
    dataset_id=z$dataset_id[1],
    C3_n=nrow(z),
    EC_identity_median=median(z$EC_identity_z,na.rm=TRUE),
    Fibroblast_median=median(z$Fibroblast_z,na.rm=TRUE),
    SMC_median=median(z$SMC_z,na.rm=TRUE),
    Pericyte_median=median(z$Pericyte_z,na.rm=TRUE),
    Myeloid_median=median(z$Myeloid_z,na.rm=TRUE),
    any_nonEC_high_pct=100*mean(z$any_nonEC_high),
    combined_nonEC_risk_pct=100*mean(z$combined_nonEC_risk)
  )
}))
write_csv(donor_summary, file.path(OUT,"05_C3_donor_level_risk.csv"))

# ------------------ CANDIDATE GENE RELATION TO RISK ------------------
candidates <- intersect(c("ANGPT2","LOX","CTHRC1","SULF1"), rownames(dat))
cand_rows <- list()
for (g in candidates) {
  x <- as.numeric(dat[g,rownames(c3md)])
  for (flag in c("Fibroblast_high","SMC_high","Pericyte_high","Myeloid_high",
                 "any_nonEC_high","combined_nonEC_risk")) {
    f <- c3md[[flag]]
    cand_rows[[length(cand_rows)+1]] <- data.frame(
      gene=g, risk_flag=flag,
      clean_n=sum(!f), risk_n=sum(f),
      clean_pct_expr=100*mean(x[!f]>0),
      risk_pct_expr=if(sum(f)>0) 100*mean(x[f]>0) else NA_real_,
      clean_mean=mean(x[!f]),
      risk_mean=if(sum(f)>0) mean(x[f]) else NA_real_,
      clean_median=median(x[!f]),
      risk_median=if(sum(f)>0) median(x[f]) else NA_real_
    )
  }
}
write_csv(do.call(rbind,cand_rows), file.path(OUT,"06_candidate_expression_risk_vs_clean_C3.csv"))

# Spearman correlations of candidate expression with program scores.
cor_rows <- list()
for (g in candidates) {
  x <- as.numeric(dat[g,rownames(c3md)])
  for (sc in c("EC_identity_z","Fibroblast_z","SMC_z","Pericyte_z","Myeloid_z")) {
    ct <- suppressWarnings(cor.test(x,c3md[[sc]],method="spearman",exact=FALSE))
    cor_rows[[length(cor_rows)+1]] <- data.frame(
      gene=g, score=sc, rho=unname(ct$estimate), PValue=ct$p.value
    )
  }
}
cor_df <- do.call(rbind,cor_rows)
cor_df$FDR <- p.adjust(cor_df$PValue,method="BH")
write_csv(cor_df,file.path(OUT,"07_candidate_score_correlations_C3.csv"))

# ------------------ DOUBLET METADATA DESCRIPTIVE ------------------
dbl_cols <- names(cmd)[grepl("doublet",names(cmd),ignore.case=TRUE)]
if (length(dbl_cols)) {
  dbl_summary <- list()
  if ("doublet_score" %in% names(cmd)) {
    for (st in CORE_STATES) {
      z <- cmd[cmd$paper_state==st,,drop=FALSE]
      dbl_summary[[length(dbl_summary)+1]] <- data.frame(
        state=st,n=nrow(z),
        median_doublet_score=median(z$doublet_score,na.rm=TRUE),
        mean_doublet_score=mean(z$doublet_score,na.rm=TRUE)
      )
    }
    write_csv(do.call(rbind,dbl_summary),file.path(OUT,"08_doublet_score_by_state_descriptive.csv"))
  }
  if ("doublet_class" %in% names(cmd)) {
    dt <- as.data.frame(table(state=cmd$paper_state,doublet_class=cmd$doublet_class),stringsAsFactors=FALSE)
    write_csv(dt,file.path(OUT,"08b_doublet_class_by_state_descriptive.csv"))
  }
}

# ------------------ DONOR-BLOCKED edgeR SENSITIVITY ------------------
make_pb <- function(counts, mdx) {
  cells <- rownames(mdx)[mdx$paper_state %in% c("C7","C3")]
  m <- mdx[cells,,drop=FALSE]
  gid <- paste(m$donor_id,m$paper_state,sep="|||")
  gf <- factor(gid)
  A <- Matrix::sparse.model.matrix(~0+gf)
  colnames(A) <- levels(gf)
  pb_counts <- counts[,cells,drop=FALSE] %*% A
  sg <- strsplit(levels(gf),"|||",fixed=TRUE)
  pm <- data.frame(
    pb_sample_id=levels(gf),
    donor_id=vapply(sg,`[[`,character(1),1),
    state=vapply(sg,`[[`,character(1),2),
    n_cells=as.integer(table(gf)),
    stringsAsFactors=FALSE
  )
  colnames(pb_counts) <- pm$pb_sample_id
  list(counts=pb_counts,meta=pm)
}

fit_edger <- function(mdx, label) {
  pb <- make_pb(counts,mdx)
  m <- pb$meta
  d7 <- m$donor_id[m$state=="C7" & m$n_cells>=5]
  d3 <- m$donor_id[m$state=="C3" & m$n_cells>=5]
  paired <- sort(intersect(d7,d3))
  sm <- m[m$donor_id %in% paired & m$state %in% c("C7","C3"),,drop=FALSE]
  sm$state <- factor(sm$state,levels=c("C7","C3"))
  sm <- sm[order(sm$donor_id,sm$state),,drop=FALSE]
  y <- pb$counts[,sm$pb_sample_id,drop=FALSE]
  if (length(paired)<3) {
    return(list(ok=FALSE, summary=data.frame(filter=label,n_paired_donors=length(paired),
                                              n_C3_cells=sum(sm$n_cells[sm$state=="C3"]),
                                              n_C7_cells=sum(sm$n_cells[sm$state=="C7"]))))
  }
  sm$condition <- factor(as.character(sm$state),levels=c("C7","C3"))
  sm$donor_factor <- factor(sm$donor_id)
  design <- model.matrix(~donor_factor+condition,data=sm)

  dge <- DGEList(counts=y)
  keepg <- filterByExpr(dge, design=design, min.count=10, min.total.count=15,
                        large.n=10, min.prop=0.7)
  dge <- dge[keepg,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,method="TMM")
  dge <- estimateDisp(dge,design,robust=TRUE)
  fit <- glmQLFit(dge,design,robust=TRUE)
  qlf <- glmQLFTest(fit,coef="conditionC3")
  tab <- topTags(qlf,n=Inf,sort.by="none")$table
  tab$gene <- rownames(tab)
  tab$FDR <- p.adjust(tab$PValue,method="BH")
  key <- merge(data.frame(gene=c("ANGPT2","LOX","CTHRC1","SULF1")),
               tab[,c("gene","logFC","logCPM","F","PValue","FDR")],
               by="gene",all.x=TRUE,sort=FALSE)
  key$filter <- label
  key$n_paired_donors <- length(paired)
  key$n_C3_cells <- sum(sm$n_cells[sm$state=="C3"])
  key$n_C7_cells <- sum(sm$n_cells[sm$state=="C7"])
  key$n_tested_genes <- nrow(tab)
  list(ok=TRUE,key=key,full=tab,meta=sm,paired=paired)
}

# Same fixed risk rule is applied to BOTH C3 and C7.
filters <- list(
  unfiltered = rep(TRUE,nrow(cmd)),
  exclude_fibroblast_high = !cmd$Fibroblast_high,
  exclude_SMC_pericyte_high = !(cmd$SMC_high | cmd$Pericyte_high),
  exclude_any_nonEC_high = !cmd$any_nonEC_high,
  exclude_combined_nonEC_risk = !cmd$combined_nonEC_risk
)

sens <- list()
fit_store <- list()
for (nm in names(filters)) {
  mdx <- cmd[filters[[nm]],,drop=FALSE]
  ans <- fit_edger(mdx,nm)
  fit_store[[nm]] <- ans
  if (isTRUE(ans$ok)) {
    sens[[nm]] <- ans$key
    write_csv(ans$meta,file.path(OUT,paste0("09_design_",nm,".csv")))
  } else {
    sens[[nm]] <- ans$summary
  }
}
sens_df <- do.call(rbind,sens)
write_csv(sens_df,file.path(OUT,"10_edgeR_keygene_contamination_sensitivity.csv"))

# HARD reproduction check for unfiltered main analysis.
# HOTFIX: 'sens' stores compact tables only; use the full fit object
# retained in fit_store for ok/paired/key fields.
u <- fit_store[["unfiltered"]]
if (!isTRUE(u$ok)) stop("Unfiltered donor-blocked model failed.")
if (length(u$paired) != 12) stop("Unfiltered model did not reproduce 12 paired donors.")

expected_key <- data.frame(
  gene=c("ANGPT2","LOX","CTHRC1","SULF1"),
  expected_logFC=c(2.982956,5.236814,3.336996,-3.245341)
)
obs <- merge(expected_key,u$key[,c("gene","logFC")],by="gene",all.x=TRUE)
obs$abs_diff <- abs(obs$logFC-obs$expected_logFC)
write_csv(obs,file.path(OUT,"11_unfiltered_keygene_reproduction_check.csv"))
if (any(!is.finite(obs$abs_diff)) || any(obs$abs_diff > 0.02)) {
  stop("Unfiltered key-gene logFC values do not reproduce frozen final analysis within tolerance.")
}

# ------------------ COMPACT REPORT ------------------
c3risk <- risk_summary[risk_summary$state=="C3",]
cat("\n===== C3 RISK SUMMARY =====\n")
print(c3risk)
cat("\n===== edgeR KEY-GENE SENSITIVITY =====\n")
print(sens_df)

report <- c(
  "CORRECTED C3 CONTAMINATION / IDENTITY AUDIT",
  "",
  paste0("Frozen object: ",EC_RDS),
  paste0("State column: ",state_col),
  "Hard QC required and passed: strict EC=6152; core=4435; C1=1440; C2=1181; C3=696; C4=624; C5=324; C7=170.",
  "",
  "Scoring:",
  "- Gene-wise z scores were computed from log-normalized RNA expression across the 4,435 strict-core EC cells.",
  "- EC/fibroblast/SMC/pericyte/myeloid program score = mean gene-wise z score for detected genes.",
  "- Non-EC-high thresholds were fixed from the C7 reference distribution at the 95th percentile.",
  "- EC-low threshold was fixed at the C7 EC-identity 25th percentile.",
  "- combined_nonEC_risk = any non-EC-high AND EC-low.",
  "- These are sensitivity flags, not proof that a cell is a contaminant.",
  "",
  "Sensitivity DE:",
  "- Risk filters were applied to BOTH C3 and C7.",
  "- Donors were reselected using >=5 cells in each state after each filter.",
  "- Raw RNA counts were re-aggregated to donor-state pseudobulk.",
  "- edgeR design: ~ donor_factor + condition (C7 reference; positive logFC = C3>C7).",
  "- filterByExpr: min.count=10, min.total.count=15, large.n=10, min.prop=0.7.",
  "",
  "Interpretation rule:",
  "- Strong robustness: ANGPT2/LOX/CTHRC1 retain positive logFC across non-EC filters and SULF1 retains negative logFC, preferably with similar effect sizes.",
  "- Concern: a key C3-up gene loses direction or most of its effect specifically after fibroblast/SMC/pericyte-risk exclusion.",
  "- FDR changes must be interpreted together with retained donor number and cell number."
)
writeLines(report,file.path(OUT,"12_METHODS_AND_INTERPRETATION.txt"))

cat("\nDONE\n")
cat("Output:",OUT,"\n")
RSCRIPT

Rscript "$R"

# Small bundle only: no RDS and no cell-by-gene matrices.
tar -czf "$BUNDLE" -C /root \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/00_HARD_QC_state_counts.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/01_marker_expression_by_state.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/02_program_scores_by_state.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/03_C7_reference_thresholds.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/04_nonEC_risk_by_state.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/05_C3_donor_level_risk.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/06_candidate_expression_risk_vs_clean_C3.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/07_candidate_score_correlations_C3.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/10_edgeR_keygene_contamination_sensitivity.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/11_unfiltered_keygene_reproduction_check.csv \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/12_METHODS_AND_INTERPRETATION.txt \
  AS_EC_C3_CONTAMINATION_AUDIT_CORRECTED_20260919/RUN_LOG.txt

echo
echo "============================================================"
echo "DONE"
echo "Upload only this small bundle:"
echo "$BUNDLE"
ls -lh "$BUNDLE"
echo "============================================================"
