# Public data acquisition pipeline

`python3 scripts/build_public_data.py --dataset GSE50084` downloads the official processed GEO
matrix, verifies a pinned SHA-256 checksum, compiles sample metadata and preserves submitted
probe-level expression. `scripts/build_bulk_cohorts.py` does the same for GSE36059 and GSE275824.

Compilation does not update diagnoses or make a cohort clinically validated. Platform-specific gene
mapping is a separate Bioconductor step. GSE50084 blood and biopsy are not assumed paired because a
safe patient linkage is unavailable in the public titles. GSE36059 nephrectomies, mixed rejection,
TCMR, ABMR and non-rejecting allografts remain distinct. GSE275824 retains public repeated-measure
patient identifiers, randomized arm and timepoint.

A checksum or sample-count mismatch stops the pipeline. Raw and derived matrices remain gitignored
and are regenerated from public sources.
