# Prespecified analysis specification

## Analysis populations

GSE50084 biopsy and blood compartments are analyzed separately. Primary discovery uses biopsies.
GSE36059 uses ABMR, TCMR and non-rejecting allografts; mixed and nephrectomy groups are excluded from
primary contrasts but retained for sensitivity analyses. GSE275824 retains every available randomized
biopsy and blocks repeated observations by public patient identifier.

## Statistical rules

- Expression: platform-specific gene mapping followed by limma moderated models.
- Multiplicity: Benjamini–Hochberg FDR within each declared contrast family.
- Enrichment: preranked Hallmark GSEA; report NES, FDR and leading edge.
- Longitudinal treatment: `(follow-up − baseline) felzartamab − (follow-up − baseline) placebo`, with
  patient blocking. Published scores use mixed-effects models with patient random intercepts.
- Replication: same-direction effect plus validation-cohort FDR; pathway replication is reported
  separately from gene replication.
- Biological replicate: sample/patient. Cells, spots and ROIs are never substituted for patients.

## Prohibited analytical shortcuts

No pooling of blood and biopsy, ABMR and TCMR, mixed and pure rejection, stable allograft and
nephrectomy, or cells across patients as independent replicates. No feature selection before data
splitting. No interpretation of GSEA as causal evidence. No clinical prediction claim from public
retrospective expression alone.
