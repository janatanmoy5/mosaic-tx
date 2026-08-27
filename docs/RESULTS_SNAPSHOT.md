# Verified public-data results snapshot

Generated locally from public GEO matrices; large derived tables are gitignored and rebuilt by scripts.

## GSE50084 discovery

- 115 samples, 29,096 probes, 19,998 Bioconductor-mapped genes.
- Biopsy DSA-positive rejection versus DSA-negative control: 1,468 genes at BH FDR < 0.05;
  CXCL11 is the leading gene and Hallmark allograft rejection is the leading pathway.
- Rejection versus DSA-positive/no-rejection: no gene survives FDR. This supports treating the
  DSA-positive/no-rejection state as a biologically active difficult control.

## GSE36059 external validation

- 411 biopsies: 65 ABMR, 35 TCMR, 3 mixed, 281 non-rejecting, and 8 nephrectomy references
  (remaining submitter categories are preserved in the compiled table).
- 21,358 mapped genes.
- ABMR versus non-rejecting: 1,340 FDR-significant genes; leading pathway IFN-gamma response.
- ABMR versus TCMR: 3,231 FDR-significant genes, supporting rejection-type specificity.

## GSE275824 randomized treatment-response validation

- 63 biopsies from 22 patients; placebo/felzartamab and baseline/week-24/week-52 retained.
- 19,434 mapped genes.
- Patient-blocked difference-in-differences: no individual gene survives FDR at week 24 or 52,
  while Hallmark pathway shifts are detected. These are pathway-level supportive findings, not
  proof of a compound-specific mechanism.

All counts must be regenerated and independently reviewed before publication. No result constitutes
clinical or experimental validation.
