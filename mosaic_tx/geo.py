"""Checksum-verified acquisition and compilation of public GEO matrices."""
from __future__ import annotations
import csv, gzip, hashlib, json, shutil, tempfile
from pathlib import Path
from urllib.request import Request, urlopen
from .io import DataContractError

def file_sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024), b""): h.update(block)
    return h.hexdigest()

def download_verified(url: str, target: Path, expected: str, force: bool=False) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        actual=file_sha256(target)
        if actual == expected: return actual
        raise DataContractError(f"cached checksum mismatch: {target}; use --force-download")
    tmp=None
    try:
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as out:
            tmp=Path(out.name)
            req=Request(url, headers={"User-Agent":"MOSAIC-Tx/0.2 research compiler"})
            with urlopen(req, timeout=120) as response: shutil.copyfileobj(response,out)
        actual=file_sha256(tmp)
        if actual != expected: raise DataContractError(f"checksum mismatch: expected {expected}, observed {actual}")
        tmp.replace(target); return actual
    finally:
        if tmp and tmp.exists(): tmp.unlink()

def _phenotype(title: str):
    rules={
      "Biopsy_No_DSA_":("biopsy","DSA_negative_control","negative","no_rejection"),
      "Blood_No_DSA_":("blood","DSA_negative_control","negative","no_rejection"),
      "Biopsy_DSA_normal_":("biopsy","DSA_positive_no_rejection","positive","no_rejection"),
      "Blood_DSA_normal_":("blood","DSA_positive_no_rejection","positive","no_rejection"),
      "Biopsy_DSA_rejection_":("biopsy","DSA_positive_rejection","positive","acute_or_chronic_rejection"),
      "Blood_DSA_rejection_":("blood","DSA_positive_rejection","positive","acute_or_chronic_rejection")}
    for prefix,value in rules.items():
        if title.startswith(prefix): return value
    raise DataContractError(f"unmapped sample title: {title}")

def compile_gse50084(source: Path, output: Path, expected_samples: int=115):
    output.mkdir(parents=True, exist_ok=True); meta={}; features=0; in_table=False; writer=None
    with gzip.open(output/"expression_probe_matrix.csv.gz","wt",newline="",encoding="utf-8") as out:
      with gzip.open(source,"rt",newline="",encoding="utf-8") as src:
       for line in src:
        line=line.rstrip("\r\n")
        if line=="!series_matrix_table_begin": in_table=True; continue
        if line=="!series_matrix_table_end": in_table=False; continue
        if in_table:
          fields=next(csv.reader([line],delimiter="\t"))
          fields=[v[1:-1] if len(v)>1 and v[0]==v[-1]=='"' else v for v in fields]
          if writer is None: writer=csv.writer(out)
          else: features+=1
          writer.writerow(fields); continue
        if line.startswith(("!Sample_title\t","!Sample_geo_accession\t","!Sample_platform_id\t")):
          fields=next(csv.reader([line],delimiter="\t")); meta[fields[0]]=fields[1:]
    titles=meta.get("!Sample_title",[]); ids=meta.get("!Sample_geo_accession",[]); platforms=meta.get("!Sample_platform_id",[])
    if not(len(titles)==len(ids)==len(platforms)==expected_samples):
      raise DataContractError(f"metadata count mismatch: {len(titles)}/{len(ids)}/{len(platforms)}, expected {expected_samples}")
    if not features: raise DataContractError("no expression features found")
    rows=[]
    for accession,title,platform in zip(ids,titles,platforms):
      compartment,group,dsa,rejection=_phenotype(title)
      rows.append({"dataset":"GSE50084","sample_accession":accession,"sample_title":title,"compartment":compartment,"analysis_group":group,"dsa_status":dsa,"rejection_status":rejection,"platform":platform,"label_review_status":"requires_clinical_review"})
    with (output/"samples.csv").open("w",newline="",encoding="utf-8") as f:
      w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    counts={}
    for row in rows:
      key=f"{row['compartment']}|{row['analysis_group']}"; counts[key]=counts.get(key,0)+1
    manifest={"schema_version":"1.0","dataset":"GSE50084","source_sha256":file_sha256(source),"sample_count":len(rows),"feature_count":features,"group_counts":counts,"expression_level":"submitted_normalized_probe_values","compiled":True,"gene_mapping_complete":False,"label_caveat":"historical GEO submitter labels"}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

def validate_compiled(output: Path, expected_samples: int=115):
    manifest=json.loads((output/"manifest.json").read_text())
    with (output/"samples.csv").open(newline="") as f: rows=list(csv.DictReader(f))
    if len(rows)!=expected_samples or manifest["sample_count"]!=expected_samples: raise DataContractError("sample count failed")
    if len({r["sample_accession"] for r in rows})!=expected_samples: raise DataContractError("duplicate accessions")
    if {r["compartment"] for r in rows}!={"biopsy","blood"}: raise DataContractError("missing compartment")
    if not (output/"expression_probe_matrix.csv.gz").exists(): raise DataContractError("missing matrix")
    return manifest
