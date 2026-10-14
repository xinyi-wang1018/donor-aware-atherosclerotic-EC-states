CURRENT MANUSCRIPT FIGURE SOURCE-DATA PACKAGE
Human atherosclerosis endothelial-cell scRNA-seq
Frozen logic: 2026-08-27 / 2026-08-31 strict-source manuscript

ABSOLUTE RULES
1. Do NOT use any old 14-donor / 165-combination result.
2. Do NOT read /root/AS_EC_MANUSCRIPT_FIGDATA_FINAL/Figure4 as manuscript evidence.
3. Strict EC discovery = GSE155512 + GSE253904 only.
4. Expected strict high-confidence EC = 6,152.
5. Expected core EC = 4,435:
   C1=1,440; C2=1,181; C3=696; C4=624; C5=324; C7=170.
6. Main C3-vs-C7 inference = donor-level paired pseudobulk, n=12:
   GSE253904=11 donors; GSE155512=1 donor.
7. Expected C3-vs-C7 edgeR:
   tested genes=4,491; FDR<0.05=1,068; C3-up (FDR<0.05 & logFC>0.5)=574.
8. Key frozen effects:
   ANGPT2 logFC=3.01; LOX=5.26; CTHRC1=3.37; SULF1=-3.21.
9. Strict arterial sensitivity:
   2,643 cells; 9 evaluable donors, all GSE253904.
   This is sensitivity analysis, NOT cross-cohort independent replication.
10. GSE159677 PA is plaque-adjacent tissue, NOT healthy normal tissue.
11. GSE247238 mapped-state order is cross-sectional state similarity, NOT lineage/time.
12. GSE218166 shows state-stage PARTIAL correspondence only:
    never label C7=Normal, C3=Early, C1/C2=Advanced.
13. All inferential dots/lines should represent donor/patient, not individual cells.
14. NanoBanana must never create statistical data. Use it only for schematic/working-model panels.

QC
- Before plotting Figure 2–4, inspect:
  EXPORT_STATUS.tsv
  Diagnostics/F2_hard_count_QC.csv
  Diagnostics/F3_donor_hard_QC.csv
  Diagnostics/F3_edgeR_count_QC.csv
  Diagnostics/F3_keygene_QC.csv
  Diagnostics/F4_arterial_*_QC.csv
- If any DO_NOT_PLOT*.txt exists in a Figure folder, stop and inspect the corresponding QC.
