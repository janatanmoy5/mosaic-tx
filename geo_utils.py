"""Standalone checksum-verified GEO acquisition and cohort compilation utilities."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path
from urllib.request import Request, urlopen


class DataContractError(ValueError):
    """Raised when a downloaded or compiled dataset violates its declared contract."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_verified(url: str, target: Path, expected: str, force: bool = False) -> str:
    """Download atomically and reject both corrupt caches and checksum mismatches."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        actual = file_sha256(target)
        if actual == expected:
            return actual
        raise DataContractError(f"cached checksum mismatch: {target}; use --force-download")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as output:
            temporary = Path(output.name)
            request = Request(url, headers={"User-Agent": "MOSAIC-Tx/0.2 research compiler"})
            with urlopen(request, timeout=120) as response:
                shutil.copyfileobj(response, output)
        actual = file_sha256(temporary)
        if actual != expected:
            raise DataContractError(f"checksum mismatch: expected {expected}, observed {actual}")
        temporary.replace(target)
        return actual
    finally:
        if temporary and temporary.exists():
            temporary.unlink()


def _phenotype(title: str) -> tuple[str, str, str, str]:
    rules = {
        "Biopsy_No_DSA_": ("biopsy", "DSA_negative_control", "negative", "no_rejection"),
        "Blood_No_DSA_": ("blood", "DSA_negative_control", "negative", "no_rejection"),
        "Biopsy_DSA_normal_": ("biopsy", "DSA_positive_no_rejection", "positive", "no_rejection"),
        "Blood_DSA_normal_": ("blood", "DSA_positive_no_rejection", "positive", "no_rejection"),
        "Biopsy_DSA_rejection_": ("biopsy", "DSA_positive_rejection", "positive", "acute_or_chronic_rejection"),
        "Blood_DSA_rejection_": ("blood", "DSA_positive_rejection", "positive", "acute_or_chronic_rejection"),
    }
    for prefix, value in rules.items():
        if title.startswith(prefix):
            return value
    raise DataContractError(f"unmapped sample title: {title}")


def compile_gse50084(source: Path, output: Path, expected_samples: int = 115) -> dict[str, object]:
    """Compile GSE50084 while preserving biopsy/blood and all three DSA states."""
    output.mkdir(parents=True, exist_ok=True)
    metadata: dict[str, list[str]] = {}
    features, in_table, writer = 0, False, None
    with gzip.open(output / "expression_probe_matrix.csv.gz", "wt", newline="", encoding="utf-8") as destination:
        with gzip.open(source, "rt", newline="", encoding="utf-8") as input_handle:
            for line in input_handle:
                line = line.rstrip("\r\n")
                if line == "!series_matrix_table_begin":
                    in_table = True
                    continue
                if line == "!series_matrix_table_end":
                    in_table = False
                    continue
                if in_table:
                    fields = next(csv.reader([line], delimiter="\t"))
                    fields = [value[1:-1] if len(value) > 1 and value[0] == value[-1] == '"' else value for value in fields]
                    if writer is None:
                        writer = csv.writer(destination)
                    else:
                        features += 1
                    writer.writerow(fields)
                elif line.startswith(("!Sample_title\t", "!Sample_geo_accession\t", "!Sample_platform_id\t")):
                    fields = next(csv.reader([line], delimiter="\t"))
                    metadata[fields[0]] = fields[1:]
    titles = metadata.get("!Sample_title", [])
    accessions = metadata.get("!Sample_geo_accession", [])
    platforms = metadata.get("!Sample_platform_id", [])
    if not (len(titles) == len(accessions) == len(platforms) == expected_samples):
        raise DataContractError("GSE50084 metadata count mismatch")
    if not features:
        raise DataContractError("no expression features found")
    records = []
    for accession, title, platform in zip(accessions, titles, platforms):
        compartment, group, dsa, rejection = _phenotype(title)
        records.append({"dataset": "GSE50084", "sample_accession": accession, "sample_title": title,
                        "compartment": compartment, "analysis_group": group, "dsa_status": dsa,
                        "rejection_status": rejection, "platform": platform,
                        "label_review_status": "requires_clinical_review"})
    with (output / "samples.csv").open("w", newline="", encoding="utf-8") as handle:
        csv_writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        csv_writer.writeheader(); csv_writer.writerows(records)
    counts: dict[str, int] = {}
    for record in records:
        key = f"{record['compartment']}|{record['analysis_group']}"
        counts[key] = counts.get(key, 0) + 1
    manifest = {"schema_version": "1.0", "dataset": "GSE50084", "source_sha256": file_sha256(source),
                "sample_count": len(records), "feature_count": features, "group_counts": counts,
                "expression_level": "submitted_normalized_probe_values", "compiled": True,
                "gene_mapping_complete": False, "label_caveat": "historical GEO submitter labels"}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def validate_compiled(output: Path, expected_samples: int = 115) -> dict[str, object]:
    manifest = json.loads((output / "manifest.json").read_text())
    with (output / "samples.csv").open(newline="") as handle:
        records = list(csv.DictReader(handle))
    if len(records) != expected_samples or manifest["sample_count"] != expected_samples:
        raise DataContractError("sample count failed")
    if len({record["sample_accession"] for record in records}) != expected_samples:
        raise DataContractError("duplicate accessions")
    if {record["compartment"] for record in records} != {"biopsy", "blood"}:
        raise DataContractError("missing compartment")
    if not (output / "expression_probe_matrix.csv.gz").exists():
        raise DataContractError("missing matrix")
    return manifest


def compile_geo(source: Path, output: Path, dataset: str, expected: int) -> dict[str, object]:
    """Compile GSE36059 or GSE275824 using cohort-specific public metadata."""
    output.mkdir(parents=True, exist_ok=True)
    metadata: dict[str, list[str]] = {}; characteristics = []
    inside, writer, features = False, None, 0
    with gzip.open(output / "expression_probe_matrix.csv.gz", "wt", newline="") as destination:
        with gzip.open(source, "rt", newline="") as input_handle:
            for line in input_handle:
                line = line.rstrip("\r\n")
                if line == "!series_matrix_table_begin": inside = True; continue
                if line == "!series_matrix_table_end": inside = False; continue
                if inside:
                    record = next(csv.reader([line], delimiter="\t"))
                    if writer is None: writer = csv.writer(destination)
                    else: features += 1
                    writer.writerow(record); continue
                if line.startswith("!Sample_"):
                    record = next(csv.reader([line], delimiter="\t")); key, values = record[0], record[1:]
                    if key == "!Sample_characteristics_ch1": characteristics.append(values)
                    elif key in ("!Sample_title", "!Sample_geo_accession", "!Sample_platform_id"): metadata[key] = values
    titles = metadata.get("!Sample_title", []); accessions = metadata.get("!Sample_geo_accession", [])
    platforms = metadata.get("!Sample_platform_id", [])
    if not (len(titles) == len(accessions) == len(platforms) == expected):
        raise DataContractError(f"{dataset} metadata count mismatch")
    characteristic_maps = [{} for _ in range(expected)]
    for values in characteristics:
        if len(values) != expected: raise DataContractError(f"{dataset} characteristic length mismatch")
        for index, value in enumerate(values):
            key, separator, content = value.partition(":")
            if separator: characteristic_maps[index][re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")] = content.strip()
    records = []
    for index, (accession, title, platform) in enumerate(zip(accessions, titles, platforms)):
        record = {"dataset": dataset, "sample_accession": accession, "sample_title": title, "platform": platform}
        record.update(characteristic_maps[index])
        if dataset == "GSE36059":
            diagnosis = next((value for key, value in characteristic_maps[index].items() if key.startswith("diagnosis_")), None)
            if diagnosis not in {"ABMR", "TCMR", "MIXED", "non-rejecting", "Nephrectomy"}:
                raise DataContractError(f"unmapped diagnosis {diagnosis}")
            record.update({"analysis_group": diagnosis, "compartment": "biopsy"})
        elif dataset == "GSE275824":
            def pick(prefix: str): return next((value for key, value in characteristic_maps[index].items() if key.startswith(prefix)), None)
            record.update({"patient_id": pick("patient"), "treatment": pick("felzartamab_placebo"),
                           "timepoint": pick("followup_baseline"), "analysis_group": pick("felzartamab_group_index")})
            if not all(record.get(key) for key in ("patient_id", "treatment", "timepoint")):
                raise DataContractError("missing trial metadata")
            record["compartment"] = "biopsy"
        records.append(record)
    fields = []
    for record in records:
        for key in record:
            if key not in fields: fields.append(key)
    with (output / "samples.csv").open("w", newline="") as handle:
        csv_writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        csv_writer.writeheader(); csv_writer.writerows(records)
    manifest = {"schema_version": "1.0", "dataset": dataset, "source_sha256": file_sha256(source),
                "sample_count": len(records), "feature_count": features, "compiled": True,
                "analysis_ready_with_original_labels": True, "label_source": "GEO submitter characteristics"}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
