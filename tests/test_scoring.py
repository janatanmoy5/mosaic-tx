import unittest

from mosaic_tx.io import DataContractError
from mosaic_tx.scoring import Weights, rank_interventions


class ScoringTests(unittest.TestCase):
    def setUp(self):
        self.disease = [
            {"patient_id": "P1", "cell_state": "endo", "feature": "IFNG", "effect": "2", "confidence": "1"},
            {"patient_id": "P1", "cell_state": "nk", "feature": "FCGR3A", "effect": "1", "confidence": "1"},
        ]
        self.safety = [
            {"intervention": "reverse", "infection_risk": "0", "nephrotoxicity": "0", "interaction_risk": "0"},
            {"intervention": "align", "infection_risk": "0", "nephrotoxicity": "0", "interaction_risk": "0"},
        ]

    def test_opposing_effects_rank_first(self):
        perturbations = [
            {"intervention": "reverse", "cell_state": "endo", "feature": "IFNG", "effect": "-2", "evidence": "1"},
            {"intervention": "reverse", "cell_state": "nk", "feature": "FCGR3A", "effect": "-1", "evidence": "1"},
            {"intervention": "align", "cell_state": "endo", "feature": "IFNG", "effect": "2", "evidence": "1"},
            {"intervention": "align", "cell_state": "nk", "feature": "FCGR3A", "effect": "1", "evidence": "1"},
        ]
        rows = rank_interventions(self.disease, perturbations, self.safety)
        self.assertEqual(rows[0]["intervention"], "reverse")
        self.assertEqual(rows[0]["reversal_score"], 1.0)
        self.assertEqual(rows[1]["reversal_score"], -1.0)

    def test_safety_penalty_changes_priority(self):
        perturbations = [
            {"intervention": name, "cell_state": state, "feature": feature, "effect": effect, "evidence": "1"}
            for name in ("reverse", "align")
            for state, feature, effect in (("endo", "IFNG", "-2"), ("nk", "FCGR3A", "-1"))
        ]
        self.safety[0]["nephrotoxicity"] = "1"
        rows = rank_interventions(self.disease, perturbations, self.safety, Weights(nephrotoxicity=2.0))
        self.assertEqual(rows[0]["intervention"], "align")

    def test_out_of_range_confidence_is_rejected(self):
        self.disease[0]["confidence"] = "1.2"
        with self.assertRaises(DataContractError):
            rank_interventions(self.disease, [], self.safety)


if __name__ == "__main__":
    unittest.main()

