import csv, gzip, tempfile, unittest
from pathlib import Path
from geo_utils import compile_gse50084, validate_compiled
class GeoTests(unittest.TestCase):
 def test_compile(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d); source=root/"matrix.gz"
   with gzip.open(source,"wt") as f:
    f.write('!Sample_title\t"Biopsy_No_DSA_1"\t"Blood_DSA_rejection_1"\n!Sample_geo_accession\t"GSM1"\t"GSM2"\n!Sample_platform_id\t"GPL6244"\t"GPL6244"\n!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\n"probe1"\t1\t2\n!series_matrix_table_end\n')
   out=root/"out"; manifest=compile_gse50084(source,out,2); validate_compiled(out,2)
   self.assertEqual(manifest["feature_count"],1)
   with (out/"samples.csv").open() as f: rows=list(csv.DictReader(f))
   self.assertEqual(rows[1]["analysis_group"],"DSA_positive_rejection")
if __name__=="__main__": unittest.main()
