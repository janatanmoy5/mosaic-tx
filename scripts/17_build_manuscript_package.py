#!/usr/bin/env python3
"""Assemble one nonredundant, auditable manuscript package from executed outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    """Read a CSV as dictionaries without changing published labels."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, records: list[dict[str, object]], fields: list[str]) -> None:
    """Write a stable-column CSV suitable for manual sample-level review."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_study_design_svg(path: Path) -> None:
    """Generate Figure 1 deterministically from the declared cohort/estimand design."""
    boxes = [
        (45, "1. DSA-state discovery", "GSE50084 | 61 biopsies", "DSA- NR | DSA+ NR | DSA+ rejection", "Estimand: rejection among DSA+ grafts"),
        (445, "2. Independent qualification", "GSE36059 | 381 biopsies", "281 NR | 35 TCMR | 65 ABMR", "Estimands: replication and specificity"),
        (845, "3. Randomized response", "GSE275824 | 22 patients", "63 biopsies | baseline, W24, W52", "Estimand: felzartamab-placebo DID"),
    ]
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1250" height="600" viewBox="0 0 1250 600">',
        '<rect width="1250" height="600" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#18313b}.title{font-size:30px;font-weight:700}.sub{font-size:16px;fill:#4c626c}.head{font-size:20px;font-weight:700}.body{font-size:16px}.gate{font-size:15px;font-weight:700;fill:#8f3c3c}</style>',
        '<text x="45" y="55" class="title">MOSAIC-Tx: prespecified cohorts, estimands, and decision gates</text>',
        '<text x="45" y="85" class="sub">All labels are public-metadata-derived; cross-platform expression is never pooled</text>',
    ]
    fills = ("#eaf5f3", "#edf3f8", "#faf3e7")
    strokes = ("#4d9c94", "#578ab5", "#b58335")
    for i, (x, heading, cohort, groups, estimand) in enumerate(boxes):
        parts.extend([
            f'<rect x="{x}" y="135" rx="12" width="345" height="285" fill="{fills[i]}" stroke="{strokes[i]}" stroke-width="2"/>',
            f'<text x="{x+25}" y="180" class="head">{heading}</text>',
            f'<text x="{x+25}" y="225" class="body">{cohort}</text>',
            f'<text x="{x+25}" y="270" class="body">{groups}</text>',
            f'<text x="{x+25}" y="330" class="body">{estimand}</text>',
            f'<text x="{x+25}" y="375" class="sub">Model: moderated effect + 95% CI + BH FDR</text>',
        ])
        if i < 2:
            parts.append(f'<path d="M {x+350} 277 L {x+390} 277" stroke="#596b72" stroke-width="4"/><path d="M {x+390} 277 l -12 -8 v 16 z" fill="#596b72"/>')
    parts.extend([
        '<rect x="45" y="465" rx="10" width="1145" height="90" fill="#f4f7f8" stroke="#bdcbd1"/>',
        '<text x="625" y="500" text-anchor="middle" class="head">Advance only if all evidence gates pass</text>',
        '<text x="625" y="532" text-anchor="middle" class="gate">DSA discrimination | external replication | ABMR-vs-TCMR specificity | randomized response | patient-level cell localization</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts) + "\n")


def sample_manifest(root: Path) -> list[dict[str, object]]:
    """Normalize every accession while retaining source labels and analysis eligibility."""
    output: list[dict[str, object]] = []
    for dataset in ("GSE50084", "GSE36059", "GSE275824"):
        for record in rows(root / "data" / "processed" / dataset / "samples.csv"):
            group = record.get("analysis_group", "")
            compartment = record.get("compartment", "")
            if dataset == "GSE50084":
                included = compartment == "biopsy"
                reason = "primary biopsy" if included else "blood excluded from biopsy analysis"
            elif dataset == "GSE36059":
                included = group in {"ABMR", "TCMR", "non-rejecting"}
                reason = "prespecified validation class" if included else f"{group} retained but excluded from primary contrast"
            else:
                included, reason = True, "randomized longitudinal biopsy"
            output.append({
                "dataset": dataset,
                "sample_accession": record.get("sample_accession", ""),
                "sample_title": record.get("sample_title", ""),
                "platform": record.get("platform", ""),
                "compartment": compartment,
                "analysis_group": group,
                "dsa_status": record.get("dsa_status", ""),
                "rejection_status": record.get("rejection_status", ""),
                "patient_id": record.get("patient_id", ""),
                "treatment": record.get("treatment", ""),
                "timepoint": record.get("timepoint", ""),
                "included_primary_analysis": str(included).upper(),
                "inclusion_reason": reason,
                "label_source": "GEO submitter metadata; see processed samples.csv",
            })
    return output


def cohort_table(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    """Count profiles and patients for each scientifically interpretable analysis stratum."""
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = {}
    for record in samples:
        if record["included_primary_analysis"] != "TRUE":
            continue
        dataset = str(record["dataset"])
        if dataset == "GSE275824":
            label = f"{record['treatment']} / {record['timepoint']}"
            role = "randomized treatment-time stratum"
        elif dataset == "GSE36059":
            label = str(record["analysis_group"])
            role = {"non-rejecting": "control", "TCMR": "disease specificity control", "ABMR": "case"}[label]
        else:
            label = str(record["analysis_group"])
            role = {"DSA_negative_control": "DSA-negative control",
                    "DSA_positive_no_rejection": "difficult DSA-positive control",
                    "DSA_positive_rejection": "DSA-positive rejection case"}[label]
        key = (dataset, label, str(record["platform"]), role)
        grouped.setdefault(key, []).append(record)
    table = []
    for (dataset, label, platform, role), group in grouped.items():
        patient_ids = {str(x["patient_id"]) for x in group if str(x["patient_id"])}
        table.append({"dataset": dataset, "analysis_stratum": label, "role": role,
                      "platform": platform, "n_profiles": len(group),
                      "n_patients": len(patient_ids) if patient_ids else "not available",
                      "analysis_unit": "patient" if dataset == "GSE275824" else "biopsy"})
    return sorted(table, key=lambda x: (x["dataset"], x["analysis_stratum"]))


def estimand_table(root: Path) -> list[dict[str, object]]:
    """Expose every primary module estimate with uncertainty and interpretation."""
    labels = {
        "nk_identity": "NK identity", "fcgr3a_adcc_effector": "FCGR3A/Fc signaling",
        "cytotoxicity": "Cytotoxicity", "ifng_chemokine": "IFN-gamma/chemokines",
        "endothelial_activation_injury": "Endothelial activation",
        "tissue_injury_repair": "Tissue injury/repair",
    }
    output = []
    for record in rows(root / "results" / "pilot" / "atlas" / "module_effects.csv"):
        if record["module"] not in labels:
            continue
        effect, statistic = float(record["logFC"]), float(record["t"])
        standard_error = abs(effect / statistic) if statistic and math.isfinite(statistic) else math.nan
        fdr = float(record["module_fdr"])
        if record["contrast"] == "ABMR_vs_TCMR" and effect < 0 and fdr < 0.05:
            interpretation = "significantly TCMR-enriched; fails ABMR-specificity gate"
        elif fdr < 0.05:
            interpretation = "supported association in declared direction" if effect > 0 else "supported reduction/opposite direction"
        else:
            interpretation = "not supported after multiplicity correction"
        output.append({
            "dataset": record["dataset"], "estimand": record["contrast"],
            "case_group": record["case_group"], "control_group": record["control_group"],
            "n_case": record["n_case"], "n_control": record["n_control"],
            "module": labels[record["module"]], "effect": f"{effect:.6g}",
            "ci95_lower": f"{effect - 1.96 * standard_error:.6g}",
            "ci95_upper": f"{effect + 1.96 * standard_error:.6g}",
            "p_value": record["P.Value"], "fdr_bh_within_contrast": record["module_fdr"],
            "interpretation": interpretation,
            "model": "limma moderated contrast; patient-blocked for GSE275824",
        })
    return output


def module_table(root: Path) -> list[dict[str, object]]:
    """Join frozen gene definitions to observed coverage on each platform."""
    config = json.loads((root / "configs" / "mechanism_modules.json").read_text())
    coverage = rows(root / "results" / "pilot" / "atlas" / "module_coverage.csv")
    return [{"dataset": x["dataset"], "module": x["module"],
             "declared_genes": ";".join(config["modules"][x["module"]]["genes"]),
             "n_declared": x["declared_genes"], "n_observed": x["observed_genes"],
             "coverage_fraction": x["coverage_fraction"],
             "interpretation": config["modules"][x["module"]]["interpretation"]} for x in coverage]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("manuscript"))
    args = parser.parse_args(); root = args.root.resolve(); out = (root / args.output).resolve()
    figures, tables = out / "figures", out / "tables"
    figures.mkdir(parents=True, exist_ok=True); tables.mkdir(parents=True, exist_ok=True)

    # Exactly three manuscript figures; older figures remain provenance artifacts outside this package.
    write_study_design_svg(figures / "Figure_1_study_design.svg")
    figure_sources = {
        "Figure_2_primary_estimands.png": root / "results" / "pilot" / "figures" / "estimand_results_main.png",
        "Figure_2_primary_estimands.pdf": root / "results" / "pilot" / "figures" / "estimand_results_main.pdf",
        "Figure_3_sample_level_verification.png": root / "results" / "pilot" / "figures" / "dataset_results_main.png",
        "Figure_3_sample_level_verification.pdf": root / "results" / "pilot" / "figures" / "dataset_results_main.pdf",
    }
    for name, source in figure_sources.items():
        if not source.is_file():
            raise FileNotFoundError(f"required executed figure is missing: {source}")
        shutil.copy2(source, figures / name)

    samples = sample_manifest(root)
    write_csv(tables / "Table_1_cohort_composition.csv", cohort_table(samples),
              ["dataset", "analysis_stratum", "role", "platform", "n_profiles", "n_patients", "analysis_unit"])
    write_csv(tables / "Table_2_primary_estimands.csv", estimand_table(root),
              ["dataset", "estimand", "case_group", "control_group", "n_case", "n_control", "module",
               "effect", "ci95_lower", "ci95_upper", "p_value", "fdr_bh_within_contrast", "interpretation", "model"])
    write_csv(tables / "Table_S1_sample_labels.csv", samples,
              ["dataset", "sample_accession", "sample_title", "platform", "compartment", "analysis_group",
               "dsa_status", "rejection_status", "patient_id", "treatment", "timepoint",
               "included_primary_analysis", "inclusion_reason", "label_source"])
    modules = module_table(root)
    write_csv(tables / "Table_S2_module_definitions_and_coverage.csv", modules,
              ["dataset", "module", "declared_genes", "n_declared", "n_observed", "coverage_fraction", "interpretation"])

    legends = """# Manuscript figure legends\n\n## Figure 1. Prespecified study design\nCohort roles, biological estimands and stop/go logic. This panel contains no results.\n\n## Figure 2. Qualification of a treatment-responsive rejection program\n(A) GSE50084 effect estimates separate DSA exposure from rejection among DSA-positive grafts. (B) Cross-cohort concordance is displayed against ABMR-versus-TCMR specificity. (C) Hallmark GSEA direction across discovery, replication, specificity and randomized treatment contrasts. (D) Patient-blocked felzartamab-minus-placebo treatment-by-time interactions for published ABMR, IFN-gamma, NK and endothelial scores. Bars are approximate 95% confidence intervals; multiplicity is controlled by BH FDR within declared families.\n\n## Figure 3. Sample-level verification\nUnsupervised PCA, individual module-score distributions and longitudinal patient trajectories. These panels verify sample structure and labeling; they are descriptive and are not substitutes for Figure 2 inferential estimates.\n"""
    (out / "FIGURE_LEGENDS.md").write_text(legends)

    artifacts = [p for p in out.rglob("*") if p.is_file()]
    manifest = {"schema_version": "1.0", "canonical_runner": "scripts/run_three_dataset_pipeline.sh",
                "figures": 3, "tables": 4,
                "artifacts": [{"path": str(p.relative_to(root)), "sha256": sha256(p)} for p in sorted(artifacts)]}
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Manuscript package created: {out} ({len(samples)} labeled public profiles)")


if __name__ == "__main__":
    main()
