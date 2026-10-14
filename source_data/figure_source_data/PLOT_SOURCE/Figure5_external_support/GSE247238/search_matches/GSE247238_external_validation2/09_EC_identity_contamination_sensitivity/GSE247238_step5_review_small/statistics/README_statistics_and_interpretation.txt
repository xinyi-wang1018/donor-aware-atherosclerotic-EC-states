GSE247238 Step 5 sensitivity edgeR completed.

Statistical unit: patient.
S1_HM_ALL and S2_HIGH require >=10 cells per patient-state; S3_STRICT_EC requires >=5 cells per patient-state.
C5 was not tested because the external dataset contained too few C5-like cells.

Contrast availability:
 sensitivity_set   contrast n_paired_patients     paired_patients
       S1_HM_ALL   C3_vs_C7                 4 CAR5;CAR6;CAR7;CAR9
       S1_HM_ALL Late_vs_C3                 4 CAR5;CAR6;CAR7;CAR9
       S1_HM_ALL Late_vs_C7                 4 CAR5;CAR6;CAR7;CAR9
         S2_HIGH   C3_vs_C7                 2           CAR5;CAR6
         S2_HIGH Late_vs_C3                 3      CAR5;CAR6;CAR7
         S2_HIGH Late_vs_C7                 3      CAR5;CAR6;CAR9
    S3_STRICT_EC   C3_vs_C7                 2           CAR6;CAR9
    S3_STRICT_EC Late_vs_C3                 3      CAR6;CAR7;CAR9
    S3_STRICT_EC Late_vs_C7                 3      CAR5;CAR6;CAR9
 minimum_cells_per_patient_state                       status
                              10                           ok
                              10                           ok
                              10                           ok
                              10 paired patients insufficient
                              10                           ok
                              10                           ok
                               5 paired patients insufficient
                               5                           ok
                               5                           ok

Core-gene direction summary:
   gene n_sets_C3_vs_C7_available n_sets_C3_vs_C7_positive
 ANGPT2                         1                        1
 CTHRC1                         1                        1
    LOX                         1                        1
  SULF1                         0                        0
 n_sets_Late_vs_C3_available n_sets_Late_vs_C3_positive
                           3                          0
                           2                          0
                           3                          0
                           3                          3
 n_sets_Late_vs_C3_negative strict_EC_C3_vs_C7_logFC strict_EC_Late_vs_C3_logFC
                          3                       NA                  -3.400873
                          2                       NA                         NA
                          3                       NA                  -5.739128
                          0                       NA                  10.028946
         direction_classification
 partial_or_set_dependent_support
 partial_or_set_dependent_support
 partial_or_set_dependent_support
  late_anchor_direction_preserved
                                                                  interpretation
 patient-level mapped-state sensitivity; not longitudinal timing or causal proof
 patient-level mapped-state sensitivity; not longitudinal timing or causal proof
 patient-level mapped-state sensitivity; not longitudinal timing or causal proof
 patient-level mapped-state sensitivity; not longitudinal timing or causal proof

Interpretation limits:
1. Direction preserved after strict filtering argues against a signal driven solely by obvious non-EC co-expression.
2. Loss of significance after filtering may reflect lower cell number as well as contamination; inspect effect direction and patient consistency.
3. This analysis does not prove lineage, longitudinal timing, or causality.
4. Non-EC risk modules are QC indicators, not direct doublet calls.
