# Prespecified analysis specification

## Analysis populations

GSE50084 biopsy and blood compartments are analyzed separately. Primary discovery uses biopsies.
GSE36059 uses ABMR, TCMR and non-rejecting allografts; mixed and nephrectomy groups are excluded from
primary contrasts but retained for sensitivity analyses. GSE275824 retains every available randomized
biopsy and blocks repeated observations by public patient identifier.

## Statistical rules

The prespecified mechanistic chain is DSA → endothelial engagement → FCGR3A-positive NK/monocyte
activation → cytotoxic and IFNG/chemokine response → microvascular endothelial injury → longitudinal
persistence or response. Operational definitions and stop rules are frozen in `MECHANISM_MAP.md` and
`configs/mechanism_modules.json`.


- Expression: platform-specific gene mapping followed by limma moderated models.
- Multiplicity: Benjamini–Hochberg FDR within each declared contrast family.
- Enrichment: preranked Hallmark GSEA; report NES, FDR and leading edge.
- Longitudinal treatment: `(follow-up − baseline) felzartamab − (follow-up − baseline) placebo`, with
  patient blocking. Published scores use mixed-effects models with patient random intercepts.
- Replication: same-direction effect plus validation-cohort FDR; pathway replication is reported
  separately from gene replication.
- Biological replicate: sample/patient. Cells, spots and ROIs are never substituted for patients.

## Arrow-level tests

- DSA is a clinical/serological exposure from metadata, never an expression-derived label.
- Antibody engagement is not directly observed; it requires DSA context plus endothelial response and is labeled indirect.
- FCGR3A signaling is scored separately in NK and monocyte patient pseudobulk.
- Cytotoxic/IFNG programs must be tested against TCMR and, where annotated, CMV/BK infection.
- Microvascular injury integrates endothelial pseudobulk with Banff MVI rather than expression alone.
- Longitudinal persistence uses treatment-by-time estimates and confidence intervals; nonsignificance alone is insufficient.
- Each arrow is classified as supported, unsupported, or not testable.

## Prohibited analytical shortcuts

No pooling of blood and biopsy, ABMR and TCMR, mixed and pure rejection, stable allograft and
nephrectomy, or cells across patients as independent replicates. No feature selection before data
splitting. No interpretation of GSEA as causal evidence. No clinical prediction claim from public
retrospective expression alone.
