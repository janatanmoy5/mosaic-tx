# Public data pipeline

Run `python scripts/build_public_data.py --dataset GSE50084` to download the official GEO processed
matrix, verify its pinned SHA-256 checksum, compile sample metadata and submitted normalized
probe-level expression, and validate the output. Downloads are atomic, cached, and gitignored.

| Compartment | DSA-negative | DSA+/no rejection | DSA+/rejection |
|---|---:|---:|---:|
| Biopsy | 20 | 13 | 28 |
| Blood | 12 | 14 | 28 |

The compiler does not reinterpret historical diagnoses, map probes to genes, combine platforms, or
perform batch correction. The rejection label combines acute and chronic rejection and must not be
called homogeneous ABMR without clinical adjudication. A clinician and curator must approve
`samples.csv`, select a pinned GPL6244 annotation, document probe aggregation, and freeze the
contrast before modeling. `model_ready` intentionally remains false until those gates pass.

GSE36059 is registered for external validation after independent phenotype mapping. GSE228300 has
only three biopsies and is registered as a mechanistic reference, not a training cohort.
