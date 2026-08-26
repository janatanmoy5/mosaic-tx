# Data contracts and governance

## Prototype tables

`disease_signature.csv`

| Column | Meaning |
|---|---|
| `patient_id` | Pseudonymous research identifier |
| `cell_state` | Versioned cell-state ontology term |
| `feature` | Versioned gene-set/pathway identifier |
| `effect` | Signed deviation from the declared reference |
| `confidence` | Reliability in [0, 1] |

`perturbations.csv` uses `intervention`, `cell_state`, `feature`, signed `effect`, and `evidence` in
[0, 1]. `safety.csv` has normalized research penalties for infection, nephrotoxicity, and
interaction risk. The prototype refuses missing required columns, empty tables, duplicate keys,
non-finite values, and out-of-range reliability scores.

## Production additions

Every observation needs: source accession; cohort and center; assay; biospecimen compartment;
processing version; reference population; dose and duration when applicable; direct/inferred status;
license; consent/data-use code; and immutable checksum. Clinical timelines require explicit units,
time zones, censoring rules, treatment start/stop times, and adjudication version.

## Governance boundary

No protected health information, controlled cohort data, credentials, model weights trained on
restricted data, or row-level clinical exports belong in Git. Use an approved compute environment,
least-privilege access, data-use review, audit logging, and disclosure control for every export.
Synthetic fixtures must use unmistakably fictional IDs and values.

