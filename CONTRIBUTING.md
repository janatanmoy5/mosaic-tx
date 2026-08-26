# Contributing

Open an issue describing the scientific question and acceptance test before a major change. New data
sources require provenance, license/data-use terms, cohort definition, leakage assessment, and a
datasheet. New model components require a baseline, ablation, grouped validation, calibration or
uncertainty analysis, and an update to the model card.

Run `python -m unittest discover -s tests -v` before submitting a pull request. Do not commit PHI,
controlled data, secrets, large assay files, generated environments, or unreviewed clinical claims.

