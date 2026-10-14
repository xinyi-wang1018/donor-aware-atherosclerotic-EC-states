options(stringsAsFactors = FALSE)

data_dir <- "/root/AS_external_raw/GSE159677"
out_file <- "/root/AS_external_results/GSE159677/GSE159677_input_audit.txt"
suffix_file <- "/root/AS_external_results/GSE159677/GSE159677_barcode_suffix_counts.csv"

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

barcode_file <- file.path(data_dir, "barcodes.tsv.gz")
feature_file <- file.path(data_dir, "features.tsv.gz")
matrix_file  <- file.path(data_dir, "matrix.mtx.gz")

required_files <- c(barcode_file, feature_file, matrix_file)

if (!all(file.exists(required_files))) {
  missing_files <- required_files[!file.exists(required_files)]
  stop(
    "缺少以下文件：\n",
    paste(missing_files, collapse = "\n")
  )
}

get_matrix_dimensions <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))

  repeat {
    line <- readLines(con, n = 1, warn = FALSE)

    if (length(line) == 0) {
      stop("无法从matrix.mtx.gz读取矩阵维度。")
    }

    if (!startsWith(line, "%")) {
      values <- as.integer(
        strsplit(trimws(line), "[[:space:]]+")[[1]]
      )

      if (length(values) < 3 || any(is.na(values[1:3]))) {
        stop("matrix.mtx.gz维度行格式异常：", line)
      }

      return(values[1:3])
    }
  }
}

barcodes <- readLines(
  gzfile(barcode_file, open = "rt"),
  warn = FALSE
)

features <- read.delim(
  gzfile(feature_file, open = "rt"),
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE
)

matrix_dims <- get_matrix_dimensions(matrix_file)

# 标准10x barcode之后是否还有样本后缀
suffix_after_10x <- sub(
  "^[ACGTN]+-[0-9]+",
  "",
  barcodes
)

suffix_after_10x[suffix_after_10x == ""] <- "<无额外后缀>"

suffix_table <- sort(
  table(suffix_after_10x),
  decreasing = TRUE
)

suffix_df <- data.frame(
  barcode_suffix = names(suffix_table),
  n_cells = as.integer(suffix_table),
  stringsAsFactors = FALSE
)

write.csv(
  suffix_df,
  suffix_file,
  row.names = FALSE
)

sink(out_file, split = TRUE)

cat("===== GSE159677 INPUT AUDIT =====\n\n")
cat("Data directory:", data_dir, "\n\n")

cat("文件大小：\n")
for (f in required_files) {
  info <- file.info(f)
  cat(
    basename(f), ":",
    round(info$size / 1024^2, 3),
    "MB\n"
  )
}

cat("\nMatrix Market维度：\n")
cat("基因/feature数量:", matrix_dims[1], "\n")
cat("细胞/barcode数量:", matrix_dims[2], "\n")
cat("非零表达值数量:", matrix_dims[3], "\n")

cat("\n配套文件行数：\n")
cat("features行数:", nrow(features), "\n")
cat("barcodes行数:", length(barcodes), "\n")
cat("features列数:", ncol(features), "\n")

cat("\n维度一致性：\n")
cat(
  "matrix行数 == features行数:",
  matrix_dims[1] == nrow(features),
  "\n"
)
cat(
  "matrix列数 == barcodes行数:",
  matrix_dims[2] == length(barcodes),
  "\n"
)

cat("\nBarcode唯一性：\n")
cat("重复barcode数:", sum(duplicated(barcodes)), "\n")
cat("唯一barcode数:", length(unique(barcodes)), "\n")

cat("\n前20个barcodes：\n")
print(head(barcodes, 20))

cat("\n后20个barcodes：\n")
print(tail(barcodes, 20))

cat("\nBarcode额外后缀及细胞数：\n")
print(head(suffix_df, 50), row.names = FALSE)

cat("\nFeatures前10行：\n")
print(head(features, 10), row.names = FALSE)

cat("\n判定提示：\n")

if (nrow(suffix_df) > 1 ||
    suffix_df$barcode_suffix[1] != "<无额外后缀>") {
  cat(
    "检测到barcode额外后缀，可能能够据此恢复样本来源。\n"
  )
} else {
  cat(
    "没有检测到barcode额外后缀，暂时无法仅凭barcode恢复6份样本来源。\n"
  )
}

cat("\n===== AUDIT COMPLETE =====\n")

sink()

cat("\n审计完成：\n")
cat(out_file, "\n")
cat(suffix_file, "\n")
