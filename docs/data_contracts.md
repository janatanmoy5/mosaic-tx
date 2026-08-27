# Data contracts

Every compiled cohort records dataset accession, sample accession, platform, original submitter title,
analysis group and source checksum. Cohort-specific fields remain separate rather than forced into a
false universal schema. GSE275824 preserves public patient, treatment and timepoint identifiers for
repeated-measure analysis.

Raw and derived matrices are gitignored and rebuilt from public sources. No protected health
information, private clinical data, credentials or restricted cohort exports belong in the repository.
