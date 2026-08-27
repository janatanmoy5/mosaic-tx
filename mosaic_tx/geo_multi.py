"""Generic GEO matrix compiler with cohort-specific, auditable phenotype maps."""
import csv,gzip,json,re
from pathlib import Path
from .geo import file_sha256
from .io import DataContractError

def compile_geo(source:Path,out:Path,dataset:str,expected:int):
 out.mkdir(parents=True,exist_ok=True); meta={}; chars=[]; inside=False; writer=None; features=0
 with gzip.open(out/'expression_probe_matrix.csv.gz','wt',newline='') as dst:
  with gzip.open(source,'rt',newline='') as src:
   for line in src:
    line=line.rstrip('\r\n')
    if line=='!series_matrix_table_begin': inside=True; continue
    if line=='!series_matrix_table_end': inside=False; continue
    if inside:
     row=next(csv.reader([line],delimiter='\t'))
     if writer is None: writer=csv.writer(dst)
     else: features+=1
     writer.writerow(row); continue
    if line.startswith('!Sample_'):
     row=next(csv.reader([line],delimiter='\t')); key=row[0]; vals=row[1:]
     if key=='!Sample_characteristics_ch1': chars.append(vals)
     elif key in ('!Sample_title','!Sample_geo_accession','!Sample_platform_id'): meta[key]=vals
 titles=meta.get('!Sample_title',[]); ids=meta.get('!Sample_geo_accession',[]); platforms=meta.get('!Sample_platform_id',[])
 if not(len(titles)==len(ids)==len(platforms)==expected): raise DataContractError(f'{dataset} metadata count mismatch')
 charmaps=[{} for _ in range(expected)]
 for vals in chars:
  if len(vals)!=expected: raise DataContractError(f'{dataset} characteristic length mismatch')
  for i,value in enumerate(vals):
   key,sep,val=value.partition(':');
   if sep: charmaps[i][re.sub(r'[^a-z0-9]+','_',key.lower()).strip('_')]=val.strip()
 rows=[]
 for i,(sid,title,platform) in enumerate(zip(ids,titles,platforms)):
  row={'dataset':dataset,'sample_accession':sid,'sample_title':title,'platform':platform}; row.update(charmaps[i])
  if dataset=='GSE36059':
   diagnosis=next((v for k,v in charmaps[i].items() if k.startswith('diagnosis_')),None)
   if diagnosis not in {'ABMR','TCMR','MIXED','non-rejecting','Nephrectomy'}: raise DataContractError(f'unmapped diagnosis {diagnosis}')
   row['analysis_group']=diagnosis; row['compartment']='biopsy'
  elif dataset=='GSE275824':
   def pick(prefix): return next((v for k,v in charmaps[i].items() if k.startswith(prefix)),None)
   row.update({'patient_id':pick('patient'),'treatment':pick('felzartamab_placebo'),'timepoint':pick('followup_baseline'),'analysis_group':pick('felzartamab_group_index')})
   if not all(row.get(k) for k in ('patient_id','treatment','timepoint')): raise DataContractError('missing trial metadata')
   row['compartment']='biopsy'
  rows.append(row)
 keys=[]
 for row in rows:
  for key in row:
   if key not in keys: keys.append(key)
 with (out/'samples.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
 manifest={'schema_version':'1.0','dataset':dataset,'source_sha256':file_sha256(source),'sample_count':len(rows),'feature_count':features,'compiled':True,'analysis_ready_with_original_labels':True,'label_source':'GEO submitter characteristics'}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n'); return manifest
