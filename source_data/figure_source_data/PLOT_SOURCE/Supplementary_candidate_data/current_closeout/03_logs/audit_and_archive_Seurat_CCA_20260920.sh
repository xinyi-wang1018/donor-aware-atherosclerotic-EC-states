#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/root"
OUT="/root/AS_EC_SEURAT_CCA_PROVENANCE_20260920"
BUNDLE="/root/AS_EC_SEURAT_CCA_PROVENANCE_20260920.tar.gz"

rm -rf "$OUT" "$BUNDLE"
mkdir -p "$OUT"/{01_candidates,02_exact_copies,03_logs,04_sessioninfo,05_report}

echo "Seurat CCA provenance audit" > "$OUT/05_report/README.txt"
echo "Goal: recover the ACTUAL successful Seurat CCA trajectory-sensitivity script/log/sessionInfo without inventing parameters." >> "$OUT/05_report/README.txt"
echo "This audit does NOT rerun integration." >> "$OUT/05_report/README.txt"

# Strong CCA evidence. We require actual Seurat integration syntax and/or CCA-labelled successful output.
PATTERN='FindIntegrationAnchors|IntegrateData|reduction[[:space:]]*=[[:space:]]*["'\'']cca["'\'']|reduction[[:space:]]*=[[:space:]]*["'\'']CCA["'\'']|Seurat.*CCA|CCA.*Slingshot|integrated.*cca|cca.*trajectory'

# Search likely source/log files.
find "$ROOT" -xdev -type f \
  \( -iname '*.R' -o -iname '*.r' -o -iname '*.Rmd' -o -iname '*.sh' -o -iname '*.txt' -o -iname '*.log' -o -iname '*.Rout' \) \
  ! -path '*/site-packages/*' ! -path '*/node_modules/*' ! -path '*/renv/*' \
  ! -path "$OUT/*" -size -20M -print0 2>/dev/null |
while IFS= read -r -d '' f; do
  if grep -Eqi "$PATTERN" "$f" 2>/dev/null; then
    printf '\n===== FILE: %s =====\n' "$f"
    grep -Ein -C 12 "$PATTERN" "$f" 2>/dev/null || true
  fi
done > "$OUT/01_candidates/CCA_exact_hits.txt"

# Candidate file inventory with a simple evidence score.
python3 - <<'PY'
from pathlib import Path
import csv, re

root = Path("/root")
out = Path("/root/AS_EC_SEURAT_CCA_PROVENANCE_20260920")
rows = []
skip_parts = ("site-packages","node_modules","renv",".cache","AS_EC_SEURAT_CCA_PROVENANCE_20260920")
exts = {".R",".r",".Rmd",".sh",".txt",".log",".Rout"}

for p in root.rglob("*"):
    try:
        if not p.is_file() or p.suffix not in exts or p.stat().st_size > 20_000_000:
            continue
        s_path = str(p)
        if any(x in s_path for x in skip_parts):
            continue
        txt = p.read_text(errors="ignore")
    except Exception:
        continue

    score = 0
    evidence = []
    tests = [
        (r"FindIntegrationAnchors", 3, "FindIntegrationAnchors"),
        (r"IntegrateData", 3, "IntegrateData"),
        (r'reduction\s*=\s*["\']cca["\']', 5, "reduction=cca"),
        (r"RunPCA", 1, "RunPCA"),
        (r"slingshot", 2, "slingshot"),
        (r"SingleCellExperiment", 1, "SingleCellExperiment"),
        (r"sessionInfo\s*\(", 1, "sessionInfo"),
        (r"\bDONE\b|completed|success", 1, "success-marker"),
        (r"future\.globals\.maxSize|RPCA", 1, "RPCA-context"),
    ]
    for pat, pts, label in tests:
        if re.search(pat, txt, flags=re.I):
            score += pts
            evidence.append(label)

    if score >= 5:
        rows.append((score, s_path, ";".join(evidence), p.stat().st_size, p.stat().st_mtime))

rows.sort(reverse=True)
with open(out/"01_candidates"/"candidate_rank.tsv","w",newline="") as f:
    w=csv.writer(f, delimiter="\t")
    w.writerow(["score","path","evidence","size_bytes","mtime_epoch"])
    w.writerows(rows)

# Copy top candidates, preserving exact bytes.
copy_dir = out/"02_exact_copies"
copy_dir.mkdir(parents=True, exist_ok=True)
import shutil, hashlib
manifest=[]
for score, s_path, evidence, size, mt in rows[:20]:
    src=Path(s_path)
    safe=s_path.replace("/root/","").replace("/","__")
    dst=copy_dir/safe
    try:
        shutil.copy2(src,dst)
        h=hashlib.sha256(dst.read_bytes()).hexdigest()
        manifest.append((score,s_path,str(dst),h,evidence))
    except Exception:
        pass

with open(out/"02_exact_copies"/"copy_manifest.tsv","w",newline="") as f:
    w=csv.writer(f, delimiter="\t")
    w.writerow(["score","source_path","copied_path","sha256","evidence"])
    w.writerows(manifest)
PY

# Search historical sessionInfo/package-version files specifically near CCA-related content.
find "$ROOT" -xdev -type f \
  \( -iname '*sessionInfo*.txt' -o -iname '*session_info*.txt' -o -iname '*sessionInfo*.Rout' -o -iname '*package*version*.txt' \) \
  ! -path "$OUT/*" -size -5M -print0 2>/dev/null |
while IFS= read -r -d '' f; do
  rel="${f#/root/}"
  safe="${rel//\//__}"
  cp -p "$f" "$OUT/04_sessioninfo/$safe" 2>/dev/null || true
done

# Current environment, explicitly NOT historical proof unless the CCA is rerun in this same invocation later.
Rscript - <<'RS' > "$OUT/04_sessioninfo/CURRENT_ENVIRONMENT_ONLY_not_historical_proof.txt" 2>&1 || true
cat("IMPORTANT: current installed environment only; do not label as historical CCA-run sessionInfo unless a run log proves it.\n\n")
print(sessionInfo())
cat("\nSelected package versions:\n")
for (p in c("Seurat","SeuratObject","harmony","slingshot","SingleCellExperiment","Matrix")) {
  v <- tryCatch(as.character(packageVersion(p)), error=function(e) "NOT_INSTALLED")
  cat(sprintf("%-25s %s\n", p, v))
}
RS

# Copy likely CCA output/log directories/files based on names.
find "$ROOT" -xdev -type f \
  \( -iname '*cca*' -o -iname '*trajectory*sens*' -o -iname '*embedding*sens*' -o -iname '*slingshot*sens*' \) \
  ! -path "$OUT/*" -size -20M -print0 2>/dev/null |
while IFS= read -r -d '' f; do
  rel="${f#/root/}"
  safe="${rel//\//__}"
  cp -p "$f" "$OUT/03_logs/$safe" 2>/dev/null || true
done

# Compact report.
{
  echo "SEURAT CCA PROVENANCE AUDIT REPORT"
  echo "=================================="
  echo
  echo "Required final evidence:"
  echo "1) exact script that successfully ran Seurat CCA integration;"
  echo "2) successful log/output proving completion;"
  echo "3) sessionInfo from that successful run, or a controlled rerun that captures sessionInfo;"
  echo "4) no substitution with RPCA (RPCA previously failed due to future.globals size limit)."
  echo
  echo "Inspect first:"
  echo "- 01_candidates/candidate_rank.tsv"
  echo "- 02_exact_copies/copy_manifest.tsv"
  echo "- 03_logs/"
  echo "- 04_sessioninfo/"
  echo
  echo "If a historical run-specific sessionInfo is absent, do NOT invent it."
  echo "Return this bundle to ChatGPT; the exact successful script can then be frozen and, if necessary, rerun once with sessionInfo capture."
} > "$OUT/05_report/CCA_PROVENANCE_REPORT.txt"

tar -czf "$BUNDLE" -C "$ROOT" "$(basename "$OUT")"

echo
echo "============================================================"
echo "DONE"
echo "Upload only:"
echo "$BUNDLE"
ls -lh "$BUNDLE"
echo "============================================================"
