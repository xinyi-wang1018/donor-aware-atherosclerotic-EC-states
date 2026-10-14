GSE247238 Step 4 patient-level mapped-state pseudobulk completed.
Primary statistical unit: patient.
Primary minimum cells per patient-state pseudobulk: 10

State edgeR contrast availability:
   contrast n_paired_patients                       status
   C3_vs_C7                 4                           ok
   C5_vs_C3                 0 paired patients insufficient
   C5_vs_C7                 0 paired patients insufficient
 Late_vs_C5                 0 paired patients insufficient
 Late_vs_C3                 4                           ok
 Late_vs_C7                 4                           ok

Candidate state-order precedence relative to SULF1:
   gene n_evaluable_patients n_earlier_than_SULF1 fraction_earlier_than_SULF1
 TGFBR3                    2                    2                   1.0000000
  SMAD1                    3                    3                   1.0000000
   SOX4                    1                    1                   1.0000000
  HIF1A                    3                    3                   1.0000000
 ANGPT2                    4                    4                   1.0000000
  NAMPT                    3                    1                   0.3333333
  THBS1                    4                    4                   1.0000000
 CTHRC1                    3                    3                   1.0000000
    LOX                    4                    4                   1.0000000
 n_same_as_SULF1 n_later_than_SULF1           state_order_support
               0                  0   limited_directional_support
               0                  0  moderate_state_order_support
               0                  0 not_supported_or_insufficient
               0                  0  moderate_state_order_support
               0                  0    strong_state_order_support
               2                  0 not_supported_or_insufficient
               0                  0    strong_state_order_support
               0                  0  moderate_state_order_support
               0                  0    strong_state_order_support
                                                                interpretation
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality
 cross-sectional mapped-state precedence; not longitudinal timing or causality

Interpretation limits:
1. Mapped-state order supports transcriptional precedence but is not true longitudinal timing.
2. A candidate earlier than SULF1 is not automatically a causal driver.
3. C4 is excluded from the arterial trajectory and retained as an independent decoy state.
4. Low-confidence cell mappings do not enter the primary pseudobulk analysis.
5. edgeR results require both effect direction and patient consistency; FDR alone is not sufficient.
