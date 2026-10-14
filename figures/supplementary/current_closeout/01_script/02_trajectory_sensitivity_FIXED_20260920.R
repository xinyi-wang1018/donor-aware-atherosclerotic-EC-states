#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
set.seed(20260920)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(SingleCellExperiment)
  library(slingshot)
})


canonical_state <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  for (st in c("C1","C2","C3","C4","C5","C7")) {
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

first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% colnames(x)]
  if (length(hit)) hit[[1]] else NA_character_
}

attach_donor <- function(md, donor_map_file) {
  dm <- read.csv(donor_map_file, check.names=FALSE, stringsAsFactors=FALSE)
  req <- c("dataset_id","sample_id","donor_id")
  if (!all(req %in% names(dm))) {
    stop("Donor map missing: ", paste(setdiff(req, names(dm)), collapse=", "))
  }
  dm <- unique(dm[, req, drop=FALSE])
  dkey <- paste(dm$dataset_id, dm$sample_id, sep="|||MAP|||")
  ckey <- paste(md$dataset_id, md$sample_id, sep="|||MAP|||")
  md$donor_id <- dm$donor_id[match(ckey, dkey)]
  if (anyNA(md$donor_id)) stop("Some cells could not be mapped to donor.")
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

write_csv <- function(x, path) {
  write.csv(as.data.frame(x), path, row.names=FALSE, quote=TRUE, na="")
}


OUT <- "/root/AS_EC_TEACHER_P0_REPAIR_20260920/02_TRAJECTORY"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

EC_RDS <- "/root/AS_EC_state_scoring/AS_EC_highconfidence_state_scored.rds"
DONOR_MAP <- "/root/AS_EC_19day_sprint/04_donor_program_validation/04_sample_to_donor_mapping.csv"

obj0 <- readRDS(EC_RDS)
md0 <- obj0@meta.data
state_col <- first_existing(md0, c("EC_cluster","EC_cluster_std","paper_state","EC_refined_state"))
if (is.na(state_col)) stop("State column not found.")

state_all <- canonical_state(md0[[state_col]])
names(state_all) <- rownames(md0)

cells <- rownames(md0)[
  md0$dataset_id %in% c("GSE155512","GSE253904") &
  state_all %in% c("C7","C3","C5","C1","C2")
]
obj0 <- subset(obj0, cells=cells)
md <- obj0@meta.data
md$paper_state <- state_all[rownames(md)]
md <- attach_donor(md, DONOR_MAP)

counts <- get_layer(obj0, "counts")
obj <- CreateSeuratObject(counts=counts, meta.data=md)
obj <- NormalizeData(obj, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=min(3000,nrow(obj)), verbose=FALSE)
obj <- ScaleData(obj, features=VariableFeatures(obj), verbose=FALSE)
obj <- RunPCA(obj, features=VariableFeatures(obj), npcs=40, verbose=FALSE)
dat <- get_layer(obj, "data")

run_traj <- function(embed, method_name) {
  cc <- rownames(embed)
  st <- as.character(obj@meta.data[cc, "paper_state"])

  dummy <- matrix(0, nrow=1, ncol=length(cc),
                  dimnames=list("dummy", cc))
  sce <- SingleCellExperiment(assays=list(dummy=dummy))
  colData(sce)$state <- factor(st)
  reducedDims(sce)$X <- as.matrix(embed)

  sce <- slingshot(
    sce,
    clusterLabels="state",
    reducedDim="X",
    start.clus="C7",
    stretch=0,
    reweight=TRUE,
    reassign=TRUE
  )

  lins <- slingLineages(sce)
  if (length(lins) == 0) stop("No lineage returned for ", method_name)

  targets <- c("C7","C3","C5","C1","C2")
  scores <- vapply(lins, function(z) sum(unique(z) %in% targets), numeric(1))
  li <- which.max(scores)

  ptmat <- slingPseudotime(sce)
  pt <- ptmat[, li]
  lineage_text <- paste(lins[[li]], collapse="->")

  med <- aggregate(pt, list(state=st), function(x) median(x, na.rm=TRUE))
  names(med)[2] <- "median_pseudotime"
  med$method <- method_name
  med$lineage <- lineage_text

  good <- is.finite(pt)
  ptn <- rep(NA_real_, length(pt))
  ptn[good] <- (pt[good] - min(pt[good])) / (max(pt[good]) - min(pt[good]) + 1e-12)

  trends <- list()
  for (g in intersect(c("ANGPT2","LOX","CTHRC1","SULF1"), rownames(dat))) {
    x <- as.numeric(dat[g, cc])
    ct <- suppressWarnings(cor.test(x[good], pt[good], method="spearman", exact=FALSE))
    bins <- cut(ptn[good], breaks=seq(0,1,length.out=11), include.lowest=TRUE, labels=FALSE)
    bm <- tapply(x[good], bins, mean, na.rm=TRUE)
    peakbin <- if (all(!is.finite(bm))) NA_integer_ else which.max(bm)

    trends[[length(trends)+1]] <- data.frame(
      method=method_name,
      gene=g,
      rho=unname(ct$estimate),
      PValue=ct$p.value,
      peak_bin=peakbin,
      peak_position=ifelse(is.na(peakbin), NA_real_, (peakbin-0.5)/10),
      stringsAsFactors=FALSE
    )
  }

  list(median=med, trend=do.call(rbind,trends), lineage=lineage_text)
}

embeds <- list()
status <- c()

# 1) Raw RNA PCA
embeds[["RNA_PCA"]] <- Embeddings(obj, "pca")[,1:20,drop=FALSE]
status <- c(status, "RNA_PCA: success")

# 2) Harmony
obj_h <- tryCatch(
  harmony::RunHarmony(
    object=obj,
    group.by.vars="dataset_id",
    reduction.use="pca",
    dims.use=1:30,
    verbose=FALSE
  ),
  error=function(e) e
)

if (!inherits(obj_h, "error")) {
  embeds[["Harmony"]] <- Embeddings(obj_h, "harmony")[,1:20,drop=FALSE]
  status <- c(status, "Harmony: success")
} else {
  status <- c(status, paste0("Harmony: FAILED - ", conditionMessage(obj_h)))
}

# 3) Seurat RPCA; fallback to CCA if RPCA fails.
split_list <- SplitObject(obj, split.by="dataset_id")
split_list <- lapply(split_list, function(x) {
  x <- NormalizeData(x, verbose=FALSE)
  x <- FindVariableFeatures(x, selection.method="vst", nfeatures=min(2500,nrow(x)), verbose=FALSE)
  x <- ScaleData(x, features=VariableFeatures(x), verbose=FALSE)
  x <- RunPCA(x, features=VariableFeatures(x), npcs=30, verbose=FALSE)
  x
})

features <- SelectIntegrationFeatures(object.list=split_list, nfeatures=2500)

run_integrated <- function(reduction_name) {
  anchors <- FindIntegrationAnchors(
    object.list=split_list,
    anchor.features=features,
    reduction=reduction_name,
    dims=1:30,
    verbose=FALSE
  )
  integ <- IntegrateData(anchorset=anchors, dims=1:30, verbose=FALSE)
  DefaultAssay(integ) <- "integrated"
  integ <- ScaleData(integ, verbose=FALSE)
  integ <- RunPCA(integ, npcs=30, verbose=FALSE)
  Embeddings(integ, "pca")[,1:20,drop=FALSE]
}

rpca <- tryCatch(run_integrated("rpca"), error=function(e) e)
if (!inherits(rpca, "error")) {
  embeds[["Seurat_RPCA"]] <- rpca
  status <- c(status, "Seurat_RPCA: success")
} else {
  status <- c(status, paste0("Seurat_RPCA: FAILED - ", conditionMessage(rpca)))
  cca <- tryCatch(run_integrated("cca"), error=function(e) e)
  if (!inherits(cca, "error")) {
    embeds[["Seurat_CCA"]] <- cca
    status <- c(status, "Seurat_CCA: success")
  } else {
    status <- c(status, paste0("Seurat_CCA: FAILED - ", conditionMessage(cca)))
  }
}

all_med <- list()
all_trend <- list()
topology <- c()

for (nm in names(embeds)) {
  ans <- tryCatch(run_traj(embeds[[nm]], nm), error=function(e) e)
  if (inherits(ans, "error")) {
    topology <- c(topology, paste0(nm, ": FAILED - ", conditionMessage(ans)))
  } else {
    all_med[[nm]] <- ans$median
    all_trend[[nm]] <- ans$trend
    topology <- c(topology, paste0(nm, ": ", ans$lineage))
  }
}

if (length(all_med) == 0) stop("All trajectory runs failed.")

meddf <- do.call(rbind, all_med)
trenddf <- do.call(rbind, all_trend)
trenddf$FDR <- p.adjust(trenddf$PValue, method="BH")

write_csv(meddf, file.path(OUT, "01_state_pseudotime_medians_by_method.csv"))
write_csv(trenddf, file.path(OUT, "02_candidate_pseudotime_trends_by_method.csv"))
writeLines(c(status, "", topology), file.path(OUT, "03_method_status_and_topology.txt"))

# State ordering based on median pseudotime
ord <- do.call(rbind, lapply(split(meddf, meddf$method), function(z) {
  z <- z[order(z$median_pseudotime),]
  data.frame(
    method=z$method[1],
    order=paste(z$state, collapse="->"),
    stringsAsFactors=FALSE
  )
}))
write_csv(ord, file.path(OUT, "04_state_ordering_by_method.csv"))

report <- c(
  "TRAJECTORY SENSITIVITY",
  "Compares Slingshot state-space ordering across raw RNA-PCA, Harmony, and an independent Seurat integration framework (RPCA; CCA fallback).",
  "C7 is retained as the same root-like reference to isolate the effect of embedding/integration framework.",
  "",
  "Method status:",
  status,
  "",
  "Topology:",
  topology,
  "",
  "Interpretation: concordant state ordering and candidate trends across embeddings support state-space robustness, but do not establish chronological disease progression."
)
writeLines(report, file.path(OUT, "05_TRAJECTORY_INTERPRETATION.txt"))

cat("Trajectory sensitivity completed.\n")
cat("Output:", OUT, "\n")
