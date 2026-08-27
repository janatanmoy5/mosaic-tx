import gzip, tempfile, unittest
from pathlib import Path
from geo_utils import compile_geo

class MultiCohortCompilerTests(unittest.TestCase):
    def _matrix(self, path, titles, accessions, characteristics):
        with gzip.open(path, "wt") as f:
            f.write("!Sample_title\t" + "\t".join(f'"{x}"' for x in titles) + "\n")
            f.write("!Sample_geo_accession\t" + "\t".join(f'"{x}"' for x in accessions) + "\n")
            f.write("!Sample_platform_id\t" + "\t".join('"GPL570"' for _ in titles) + "\n")
            for key, values in characteristics:
                f.write("!Sample_characteristics_ch1\t" + "\t".join(f'"{key}: {x}"' for x in values) + "\n")
            f.write('!series_matrix_table_begin\n"ID_REF"\t' + "\t".join(f'"{x}"' for x in accessions) + '\n"p1"\t' + "\t".join("1" for _ in titles) + "\n!series_matrix_table_end\n")
    def test_gse36059_keeps_nephrectomy_separate(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); src=root/"x.gz"
            self._matrix(src,["ABMR1","Nephrectomy1"],["G1","G2"],[("diagnosis (tcmr,abmr,mixed,non-rejecting)",["ABMR","Nephrectomy"])])
            result=compile_geo(src,root/"out","GSE36059",2)
            self.assertEqual(result["sample_count"],2)
    def test_treatment_requires_patient_arm_and_time(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); src=root/"x.gz"; vals=["1","2"]
            self._matrix(src,["f1","f2"],["G1","G2"],[("patient",vals),("felzartamab (placebo, felzartamab)",["Placebo","Felzartamab"]),("followup (baseline, week24, week52)",["Baseline","Week24"]),("felzartamab_group (index_placebo, fu1_felzartamab)",["Index_Placebo","FU1_Felzartamab"])])
            result=compile_geo(src,root/"out","GSE275824",2)
            self.assertEqual(result["sample_count"],2)
if __name__ == "__main__": unittest.main()
