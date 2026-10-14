library(Seurat)
library(harmony)
library(future)

plan(multicore, workers = 32)
options(future.globals.maxSize = 300 * 1024^3)

cat("Step 1: Load SCT objects\n")
obj_list <- readRDS("/root/all_53_sct_objs.rds")

cat("Object number:", length(obj_list), "\n")
cat("Total cells:", sum(sapply(obj_list, ncol)), "\n")

for (i in seq_along(obj_list)) {
  DefaultAssay(obj_list[[i]]) <- "SCT"
}

cat("Step 2: Manually select integration features\n")
features <- names(
  sort(
    table(unlist(lapply(obj_list, VariableFeatures))),
    decreasing = TRUE
  )
)[1:3000]

features <- as.character(features)
cat("Feature number:", length(features), "\n")

cat("Step 3: Prep SCT integration\n")
obj_list <- PrepSCTIntegration(
  object.list = obj_list,
  anchor.features = features
)

cat("Step 4: Run PCA for each sample\n")
obj_list <- lapply(obj_list, function(x) {
  RunPCA(
    x,
    features = features,
    npcs = 50,
    verbose = FALSE
  )
})

cat("Step 5: Find integration anchors using RPCA\n")
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:50
)

cat("Step 6: Integrate data\n")
combined <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:50
)

saveRDS(combined, "/root/AS_53samples_integrated_SCT_only.rds")

cat("Step 7: PCA\n")
DefaultAssay(combined) <- "integrated"

combined <- RunPCA(
  combined,
  npcs = 50,
  verbose = FALSE
)

cat("Step 8: Harmony\n")
meta_cols <- colnames(combined@meta.data)

batch_col <- if ("dataset_id" %in% meta_cols) {
  "dataset_id"
} else if ("dataset" %in% meta_cols) {
  "dataset"
} else {
  "orig.ident"
}

cat("Using batch column:", batch_col, "\n")

combined <- RunHarmony(
  combined,
  group.by.vars = batch_col,
  reduction = "pca",
  assay.use = "integrated"
)

cat("Step 9: UMAP and clustering\n")
combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:30
)

combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:30
)

combined <- FindClusters(
  combined,
  resolution = 0.5
)

cat("Step 10: Save final object\n")
saveRDS(
  combined,
  "/root/AS_53samples_integrated_harmony_umap_clustered.rds"
)

cat("ALL DONE\n")
