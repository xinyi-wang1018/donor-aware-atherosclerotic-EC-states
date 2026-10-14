#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# Collect source materials needed to write the Methods section
# Output:
#   /root/AS_methods_source_bundle_20260805.tar.gz
# The script does NOT copy large expression matrices or RDS objects.
# It only copies relevant scripts, small metadata/result tables, logs,
# and generates compact summaries of selected Seurat objects.
# ============================================================

OUT="/root/AS_methods_source_bundle_20260805"
ARCHIVE="/root/AS_methods_source_bundle_20260805.tar.gz"

rm -rf "$OUT" "$ARCHIVE"
mkdir -p \
  "$OUT/01_scripts" \
  "$OUT/02_small_tables" \
  "$OUT/03_environment" \
  "$OUT/04_object_summary" \
  "$OUT/05_manifest"

LOG="$OUT/03_environment/collection.log"
exec > >(tee -a "$LOG") 2>&1

echo "Collection started: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"

MANIFEST="$OUT/05_manifest/copied_files.tsv"
printf "category\toriginal_path\tcopied_path\tsize_bytes\n" > "$MANIFEST"
SKIPPED="$OUT/05_manifest/skipped_possible_secret_files.txt"
: > "$SKIPPED"

# Paths commonly containing software environments or caches.
EXCLUDE_PATHS=(
  "/root/miniconda3/"
  "/root/anaconda3/"
  "/root/mambaforge/"
  "/root/.cache/"
  "/root/.local/lib/"
  "/root/R/"
  "/root/.conda/"
  "/root/.vscode-server/"
  "/root/.rstudio/"
  "/root/.ssh/"
)

is_excluded_path() {
  local src="$1"
  for p in "${EXCLUDE_PATHS[@]}"; do
    [[ "$src" == "$p"* ]] && return 0
  done
  return 1
}

contains_possible_secret() {
  local src="$1"
  grep -Eqi \
    'OPENGWAS_JWT|OPENAI_API_KEY|API[_-]?KEY[[:space:]]*[:=]|PASSWORD[[:space:]]*[:=]|Authorization:[[:space:]]*Bearer|SECRET[[:space:]]*[:=]' \
    "$src" 2>/dev/null
}

copy_with_path() {
  local category="$1"
  local src="$2"
  local base="$OUT/$category"
  local rel="${src#/}"
  local dst="$base/$rel"

  is_excluded_path "$src" && return 0

  if contains_possible_secret "$src"; then
    echo "$src" >> "$SKIPPED"
    return 0
  fi

  mkdir -p "$(dirname "$dst")"
  cp -p "$src" "$dst"
  printf "%s\t%s\t%s\t%s\n" \
    "$category" "$src" "$dst" "$(stat -c '%s' "$src" 2>/dev/null || echo NA)" \
    >> "$MANIFEST"
}

echo
echo "1. Recording software environment"

{
  echo "===== OS ====="
  uname -a || true
  echo
  echo "===== R ====="
  R --version 2>&1 || true
  echo
  echo "===== Python ====="
  python3 --version 2>&1 || true
} > "$OUT/03_environment/system_versions.txt"

cat > "$OUT/03_environment/package_versions.R" <<'RSCRIPT'
pkgs <- c(
  "Seurat", "SeuratObject", "harmony", "scDblFinder",
  "SingleCellExperiment", "slingshot", "tradeSeq",
  "edgeR", "limma", "DESeq2",
  "clusterProfiler", "org.Hs.eg.db", "fgsea",
  "liana", "nichenetr", "CellChat",
  "Matrix", "dplyr", "data.table", "ggplot2"
)

con <- file("r_package_versions.txt", open = "wt")
sink(con)
cat("R session information\n")
cat("=====================\n")
print(sessionInfo())
cat("\nSelected package versions\n")
cat("=========================\n")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT_INSTALLED")
  cat(sprintf("%-25s %s\n", p, v))
}
sink()
close(con)
RSCRIPT

(
  cd "$OUT/03_environment"
  Rscript package_versions.R || true
)

{
  python3 - <<'PY'
import importlib
import platform

mods = [
    "scanpy", "anndata", "numpy", "pandas", "scipy",
    "sklearn", "harmonypy", "scvi", "matplotlib"
]

print("Python:", platform.python_version())
for mod in mods:
    try:
        m = importlib.import_module(mod)
        print(f"{mod:20s} {getattr(m, '__version__', 'VERSION_UNKNOWN')}")
    except Exception:
        print(f"{mod:20s} NOT_INSTALLED")
PY
} > "$OUT/03_environment/python_package_versions.txt" 2>&1 || true

echo
echo "2. Searching for analysis scripts"

SCRIPT_NAME_REGEX='AS|ather|GSE155512|GSE179159|GSE224273|GSE234077|GSE253904|GSE159677|GSE247238|endothelial|(^|[_-])EC([_-]|$)|Seurat|Harmony|scDblFinder|Slingshot|pseudotime|trajectory|pseudobulk|LIANA|NicheNet|annotation|candidate|evidence|state|mapping|contamination|QC|doublet|Step[0-9]'
SCRIPT_CONTENT_REGEX='CreateSeuratObject|NormalizeData|SCTransform|FindVariableFeatures|RunPCA|RunHarmony|HarmonyMatrix|FindNeighbors|FindClusters|RunUMAP|FindAllMarkers|FindMarkers|AddModuleScore|scDblFinder|slingshot|SlingshotDataSet|tradeSeq|glmQLFit|glmQLFTest|edgeR|liana|nichenet|FindTransferAnchors|MapQuery|TransferData|pseudobulk'

while IFS= read -r -d '' f; do
  is_excluded_path "$f" && continue
  if [[ "$f" =~ $SCRIPT_NAME_REGEX ]] || grep -Eqi "$SCRIPT_CONTENT_REGEX" "$f" 2>/dev/null; then
    copy_with_path "01_scripts" "$f"
  fi
done < <(
  find /root -xdev -type f \
    \( -iname "*.R" -o -iname "*.Rmd" -o -iname "*.py" -o -iname "*.sh" -o -iname "*.ipynb" \) \
    -size -15M -print0 2>/dev/null
)

echo
echo "3. Collecting small metadata, QC, annotation and validation tables"

TABLE_REGEX='metadata|sample|patient|clinical|phenotype|group|design|mapping|barcode|QC|quality|doublet|retain|cluster|annotation|marker|module|score|pseudobulk|state|trajectory|slingshot|candidate|evidence|summary|contrast|liana|nichenet|ligand|receptor|contamination|identity|sensitivity'

while IFS= read -r -d '' f; do
  is_excluded_path "$f" && continue
  if [[ "$f" =~ $TABLE_REGEX ]]; then
    copy_with_path "02_small_tables" "$f"
  fi
done < <(
  find /root -xdev -type f \
    \( -iname "*.csv" -o -iname "*.tsv" -o -iname "*.txt" -o -iname "*.xlsx" -o -iname "*.xls" \) \
    -size -25M -print0 2>/dev/null
)

echo
echo "4. Extracting key method lines from collected scripts"

KEY_LINES="$OUT/05_manifest/key_method_lines.txt"
: > "$KEY_LINES"

if command -v rg >/dev/null 2>&1; then
  rg -n --no-heading \
    'CreateSeuratObject|nFeature_RNA|nCount_RNA|percent\.mt|PercentageFeatureSet|subset\(|scDblFinder|NormalizeData|SCTransform|FindVariableFeatures|ScaleData|RunPCA|RunHarmony|FindNeighbors|FindClusters|RunUMAP|FindAllMarkers|FindMarkers|AddModuleScore|slingshot|tradeSeq|glmQLFit|glmQLFTest|filterByExpr|calcNormFactors|liana|nichenet|FindTransferAnchors|MapQuery|TransferData|BH|p\.adjust|FDR|logFC|log2FC' \
    "$OUT/01_scripts" > "$KEY_LINES" || true
else
  grep -RniE \
    'CreateSeuratObject|nFeature_RNA|nCount_RNA|percent\.mt|PercentageFeatureSet|subset\(|scDblFinder|NormalizeData|SCTransform|FindVariableFeatures|ScaleData|RunPCA|RunHarmony|FindNeighbors|FindClusters|RunUMAP|FindAllMarkers|FindMarkers|AddModuleScore|slingshot|tradeSeq|glmQLFit|glmQLFTest|filterByExpr|calcNormFactors|liana|nichenet|FindTransferAnchors|MapQuery|TransferData|BH|p\.adjust|FDR|logFC|log2FC' \
    "$OUT/01_scripts" > "$KEY_LINES" || true
fi

echo
echo "5. Locating final Seurat/RDS objects for compact summaries"

OBJECT_PATHS="$OUT/04_object_summary/object_paths.txt"
: > "$OBJECT_PATHS"

MAIN_RDS="/root/AS_53samples_geneidfixed_RNA_harmony_cleanADT/AS_53samples_geneidfixed_cleanADT_RNA_harmony_umap_clustered.rds"
[[ -f "$MAIN_RDS" ]] && echo "$MAIN_RDS" >> "$OBJECT_PATHS"

find /root -xdev -type f -iname "*.rds" \
  \( -path "*AS_major_annotation_current33*" -o \
     -iname "*final_manual_reviewed.rds" -o \
     -iname "*EC*final*.rds" -o \
     -iname "*EC*annot*.rds" -o \
     -iname "*slingshot*.rds" -o \
     -iname "*pseudotime*.rds" \) \
  2>/dev/null | sort -u | head -n 5 >> "$OBJECT_PATHS" || true

sort -u "$OBJECT_PATHS" -o "$OBJECT_PATHS"

cat > "$OUT/04_object_summary/summarize_objects.R" <<'RSCRIPT'
suppressPackageStartupMessages({
  ok_seurat <- requireNamespace("Seurat", quietly = TRUE)
  ok_seurat_object <- requireNamespace("SeuratObject", quietly = TRUE)
})

paths <- readLines("object_paths.txt", warn = FALSE)
paths <- unique(paths[file.exists(paths)])

safe_write <- function(x, file) {
  writeLines(capture.output(x), con = file)
}

summarize_one <- function(path, idx) {
  out_txt <- sprintf("object_%02d_summary.txt", idx)
  con <- file(out_txt, open = "wt")
  sink(con)
  on.exit({
    sink()
    close(con)
  }, add = TRUE)

  cat("PATH:", path, "\n")
  cat("FILE_SIZE_BYTES:", file.info(path)$size, "\n\n")

  obj <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(obj, "error")) {
    cat("READ_ERROR:", conditionMessage(obj), "\n")
    return(invisible(NULL))
  }

  cat("CLASS:\n")
  print(class(obj))

  if (inherits(obj, "Seurat")) {
    cat("\nDIMENSIONS (features x cells):\n")
    print(dim(obj))

    cat("\nASSAYS:\n")
    print(names(obj@assays))

    cat("\nDEFAULT_ASSAY:\n")
    print(tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NA))

    cat("\nREDUCTIONS:\n")
    print(names(obj@reductions))

    cat("\nGRAPHS:\n")
    print(names(obj@graphs))

    cat("\nIDENTITY COUNTS:\n")
    print(sort(table(SeuratObject::Idents(obj)), decreasing = TRUE))

    meta <- obj@meta.data
    cat("\nMETADATA COLUMNS:\n")
    print(colnames(meta))

    candidate_cols <- c(
      "dataset_id", "dataset", "orig.ident",
      "sample_id", "sample", "patient_id", "patient",
      "condition", "disease", "group", "region",
      "seurat_clusters", "celltype", "cell_type",
      "major_celltype", "major_annotation", "EC_state", "state"
    )

    present <- intersect(candidate_cols, colnames(meta))
    for (nm in present) {
      cat("\nCOUNTS FOR:", nm, "\n")
      print(sort(table(meta[[nm]], useNA = "ifany"), decreasing = TRUE))
    }

    if ("commands" %in% slotNames(obj)) {
      cat("\nSEURAT COMMAND LOG:\n")
      for (nm in names(obj@commands)) {
        cmd <- obj@commands[[nm]]
        cat("\n---", nm, "---\n")
        cat("TIME:", tryCatch(as.character(cmd@time.stamp), error = function(e) NA), "\n")
        cat("ASSAY USED:", tryCatch(as.character(cmd@assay.used), error = function(e) NA), "\n")
        cat("PARAMETERS:\n")
        print(tryCatch(cmd@params, error = function(e) "PARAMETERS_UNAVAILABLE"))
      }
    }

    # Compact metadata counts as CSV files.
    for (nm in present) {
      tab <- as.data.frame(table(meta[[nm]], useNA = "ifany"))
      colnames(tab) <- c(nm, "n_cells")
      write.csv(
        tab,
        sprintf("object_%02d_%s_counts.csv", idx, gsub("[^A-Za-z0-9_]+", "_", nm)),
        row.names = FALSE
      )
    }
  } else {
    cat("\nSTRUCTURE:\n")
    str(obj, max.level = 2)
  }

  rm(obj)
  invisible(gc())
}

for (i in seq_along(paths)) {
  tryCatch(
    summarize_one(paths[[i]], i),
    error = function(e) {
      writeLines(
        paste("SUMMARY_ERROR:", conditionMessage(e)),
        sprintf("object_%02d_error.txt", i)
      )
    }
  )
}
RSCRIPT

(
  cd "$OUT/04_object_summary"
  Rscript summarize_objects.R || true
)

echo
echo "6. Searching for possible residual secrets in copied materials"

SECRET_REPORT="$OUT/05_manifest/POSSIBLE_SECRETS_REVIEW_BEFORE_UPLOAD.txt"
grep -RniE \
  'OPENGWAS_JWT|OPENAI_API_KEY|API[_-]?KEY|PASSWORD|Authorization:[[:space:]]*Bearer|SECRET[[:space:]]*[:=]' \
  "$OUT/01_scripts" "$OUT/02_small_tables" \
  > "$SECRET_REPORT" 2>/dev/null || true

cat > "$OUT/README_UPLOAD.txt" <<'EOF'
This bundle is intended only for drafting the Methods section.

Before uploading:
1. Open 05_manifest/POSSIBLE_SECRETS_REVIEW_BEFORE_UPLOAD.txt.
2. If it contains any real token, password, JWT or API key, delete the affected copied file.
3. Do not upload .Renviron, .Rhistory, .bash_history, SSH keys or database credentials.
4. The bundle does not include large RDS objects or raw count matrices.

Primary file to upload:
  /root/AS_methods_source_bundle_20260805.tar.gz
EOF

echo
echo "7. Creating final archive"

tar -czf "$ARCHIVE" -C /root "$(basename "$OUT")"

echo
echo "Collection completed: $(date)"
echo "Please review this file first:"
echo "  $OUT/05_manifest/POSSIBLE_SECRETS_REVIEW_BEFORE_UPLOAD.txt"
echo
echo "Then download and upload:"
echo "  $ARCHIVE"
ls -lh "$ARCHIVE"
