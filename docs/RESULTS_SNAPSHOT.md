# Verified public-data results snapshot

All values below were regenerated locally from public GEO files. Derived matrices are gitignored.

## GSE50084: DSA-state discovery

- 115 total profiles; primary biopsy groups: 20 DSA-negative/no-rejection, 13 DSA-positive/no-rejection,
  and 28 DSA-positive/rejection.
- 29,096 probes; 19,998 genes after Bioconductor mapping and probe averaging.
- Rejection versus DSA-negative: 1,468 genes at BH FDR < 0.05; CXCL11 is the leading gene and
  Hallmark allograft rejection the leading pathway.
- Rejection versus DSA-positive/no-rejection: no gene passes FDR. This null result argues against
  treating DSA-positive/no-rejection tissue as an inert control.

## GSE36059: external replication and specificity

- 411 biopsies compiled: 65 ABMR, 35 TCMR, 281 non-rejecting, 22 mixed and 8 nephrectomy references,
  with original additional metadata preserved.
- 21,358 mapped genes.
- ABMR versus non-rejecting: 1,340 FDR-significant genes; IFN-gamma response leads pathway results.
- ABMR versus TCMR: 3,231 FDR-significant genes; CCL18 leads the ranked table. This supports
  rejection-type differences but does not prove ABMR-specific causality.

## GSE275824: randomized longitudinal response

- 63 biopsies from 22 patients; placebo/felzartamab and baseline/week-24/week-52 retained.
- Genome-wide patient-blocked difference-in-differences: no single gene survives FDR at week 24 or
  week 52, while pathway-level shifts are detected.
- Mixed-effects analysis of six study-provided molecular scores found four significant week-24
  treatment interactions after BH correction: ABMR activity (estimate −0.798, FDR 0.0050), ABMR
  IFN-gamma (−0.866, FDR 0.0135), ABMR NK (−0.522, FDR 0.0193), and NK burden (−0.579, FDR 0.0341).
- ABMR endothelial and endothelial-associated scores showed no week-24 interaction. No declared
  score had a significant week-52 interaction.

This pattern motivates—but does not prove—the hypothesis of transiently suppressible immune activity
with less responsive endothelial injury. Independent review and single-cell localization remain required.
