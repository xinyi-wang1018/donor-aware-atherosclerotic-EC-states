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
  if (!all(req %in% names(dm))) stop("Donor map missing required columns: ", paste(setdiff(req,names(dm)), collapse=", "))
  dm <- unique(dm[,req,drop=FALSE])
  dkey <- paste(dm$dataset_id, dm$sample_id, sep="|||MAP|||")
  ckey <- paste(md$dataset_id, md$sample_id, sep="|||MAP|||")
  md$donor_id <- dm$donor_id[match(ckey,dkey)]
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

write_csv <- function(x, path) write.csv(as.data.frame(x), path, row.names=FALSE, quote=TRUE, na="")


OUT <- "/root/AS_EC_TEACHER_FINAL_P0_20260920/02_TRAJECTORY"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
EC_RDS <- "/root/AS_EC_state_scoring/AS_EC_highconfidence_state_scored.rds"
DONOR_MAP <- "/root/AS_EC_19day_sprint/04_donor_program_validation/04_sample_to_donor_mapping.csv"

obj0 <- readRDS(EC_RDS)
md0 <- obj0@meta.data
state_col <- first_existing(md0,c("EC_cluster","EC_cluster_std","paper_state","EC_refined_state"))
if (is.na(state_col)) stop("State column not found.")
state_all <- canonical_state(md0[[state_col]]); names(state_all) <- rownames(md0)
cells <- rownames(md0)[md0$dataset_id %in% c("GSE155512","GSE253904") & state_all %in% c("C7","C3","C5","C1","C2")]
obj0 <- subset(obj0,cells=cells)
md <- obj0@meta.data
md$paper_state <- state_all[rownames(md)]
md <- attach_donor(md,DONOR_MAP)
counts <- get_layer(obj0,"counts")

obj <- CreateSeuratObject(counts=counts,meta.data=md)
obj <- NormalizeData(obj,verbose=FALSE)
obj <- FindVariableFeatures(obj,nfeatures=min(3000,nrow(obj)),verbose=FALSE)
obj <- ScaleData(obj,features=VariableFeatures(obj),verbose=FALSE)
obj <- RunPCA(obj,features=VariableFeatures(obj),npcs=40,verbose=FALSE)

run_traj <- function(embed, cells, method, dat) {
  st <- obj@meta.data[cells,"paper_state"]
  sce <- SingleCellExperiment(assays=list(dummy=matrix(0,nrow=1,ncol=length(cells))))
  colnames(sce) <- cells
  colData(sce)$state <- factor(st)
  reducedDims(sce)$X <- as.matrix(embed[cells,,drop=FALSE])
  sce <- slingshot(sce, clusterLabels="state", reducedDim="X",
                   start.clus="C7", stretch=0, reweight=TRUE, reassign=TRUE)
  lins <- slingLineages(sce)
  targets <- c("C7","C3","C5","C1","C2")
  scores <- vapply(lins,function(z) sum(unique(z) %in% targets),numeric(1))
  li <- which.max(scores)
  pt <- slingPseudotime(sce)[,li]
  top <- paste(lins[[li]],collapse="->")
  med <- aggregate(pt, list(state=st), function(x) median(x,na.rm=TRUE))
  names(med)[2] <- "median_pseudotime"
  med$method <- method
  med$lineage <- top

  trend <- list()
  good <- is.finite(pt)
  ptn <- (pt-min(pt[good]))/(max(pt[good])-min(pt[good])+1e-12)
  for (g in intersect(c("ANGPT2","LOX","CTHRC1","SULF1"),rownames(dat))) {
    x <- as.numeric(dat[g,cells])
    ct <- suppressWarnings(cor.test(x[good],pt[good],method="spearman",exact=FALSE))
    bins <- cut(ptn[good], breaks=seq(0,1,length.out=11), include.lowest=TRUE, labels=FALSE)
    bm <- tapply(x[good],bins,mean,na.rm=TRUE)
    peakbin <- if(all(!is.finite(bm))) NA_integer_ else which.max(bm)
    trend[[length(trend)+1]] <- data.frame(method=method,gene=g,rho=unname(ct$estimate),
      PValue=ct$p.value,peak_bin=peakbin,peak_position=ifelse(is.na(peakbin),NA,(peakbin-0.5)/10))
  }
  list(median=med,trend=do.call(rbind,trend),lineage=top)
}

dat <- get_layer(obj,"data")
embeds <- list()
embeds[["RNA_PCA"]] <- Embeddings(obj,"pca")[,1:20,drop=FALSE]

obj_h <- RunHarmony(obj, group.by.vars="dataset_id", reduction="pca", dims.use=1:30,
                    assay.use="RNA", verbose=FALSE)
embeds[["Harmony"]] <- Embeddings(obj_h,"harmony")[,1:20,drop=FALSE]

rpca_status <- "not_run"
rpca_embed <- NULL
try({
  spl <- SplitObject(obj,split.by="dataset_id")
  spl <- lapply(spl,function(x){
    x <- NormalizeData(x,verbose=FALSE)
    x <- FindVariableFeatures(x,nfeatures=min(2500,nrow(x)),verbose=FALSE)
    x <- ScaleData(x,features=VariableFeatures(x),verbose=FALSE)
    x <- RunPCA(x,features=VariableFeatures(x),npcs=30,verbose=FALSE)
    x
  })
  feats <- SelectIntegrationFeatures(object.list=spl,nfeatures=2500)
  anchors <- FindIntegrationAnchors(object.list=spl,anchor.features=feats,reduction="rpca",dims=1:30,verbose=FALSE)
  integ <- IntegrateData(anchorset=anchors,dims=1:30,verbose=FALSE)
  DefaultAssay(integ) <- "integrated"
  integ <- ScaleData(integ,verbose=FALSE)
  integ <- RunPCA(integ,npcs=30,verbose=FALSE)
  rpca_embed <- Embeddings(integ,"pca")[,1:20,drop=FALSE]
  rpca_status <- "success"
}, silent=TRUE)
if (!is.null(rpca_embed)) embeds[["Seurat_RPCA"]] <- rpca_embed

allmed <- list(); alltrend <- list(); topology <- c()
for (nm in names(embeds)) {
  ans <- tryCatch(run_traj(embeds[[nm]], rownames(embeds[[nm]]), nm, dat),
                  error=function(e) e)
  if (inherits(ans,"error")) {
    topology <- c(topology,paste(nm,"FAILED:",conditionMessage(ans)))
  } else {
    allmed[[nm]] <- ans$median
    alltrend[[nm]] <- ans$trend
    topology <- c(topology,paste(nm,":",ans$lineage))
  }
}
meddf <- do.call(rbind,allmed)
trenddf <- do.call(rbind,alltrend)
if (!is.null(trenddf) && nrow(trenddf)) trenddf$FDR <- p.adjust(trenddf$PValue,method="BH")
write_csv(meddf,file.path(OUT,"01_state_pseudotime_medians_by_method.csv"))
write_csv(trenddf,file.path(OUT,"02_candidate_pseudotime_trends_by_method.csv"))
writeLines(c(paste0("RPCA status: ",rpca_status),"",topology),file.path(OUT,"03_topology_by_method.txt"))

wide <- reshape(meddf[,c("state","method","median_pseudotime")],idvar="state",timevar="method",direction="wide")
write_csv(wide,file.path(OUT,"04_state_ordering_comparison.csv"))
cat("Trajectory sensitivity done:",OUT,"\n")
