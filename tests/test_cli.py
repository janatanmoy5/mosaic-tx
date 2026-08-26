import json
from pathlib import Path
import tempfile
import unittest

from mosaic_tx.cli import main


class CliTests(unittest.TestCase):
    def test_demo_run_emits_rankings_and_audit(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "rankings.csv"
            audit = Path(directory) / "audit.json"
            code = main([
                "rank",
                "--disease", str(root / "examples/synthetic/disease_signature.csv"),
                "--perturbations", str(root / "examples/synthetic/perturbations.csv"),
                "--safety", str(root / "examples/synthetic/safety.csv"),
                "--output", str(output),
                "--audit", str(audit),
            ])
            self.assertEqual(code, 0)
            self.assertIn("priority_score", output.read_text())
            self.assertTrue(json.loads(audit.read_text())["research_use_only"])


if __name__ == "__main__":
    unittest.main()

