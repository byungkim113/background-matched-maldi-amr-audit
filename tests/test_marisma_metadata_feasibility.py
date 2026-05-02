import unittest

from scripts.marisma_metadata_feasibility_audit import markdown_table


class MarismaMetadataFeasibilityTests(unittest.TestCase):
    def test_markdown_table_does_not_require_tabulate(self):
        rows = [
            {"drug": "Ciprofloxacin", "n_s": 10, "n_r": 5},
            {"drug": "Cefotaxime", "n_s": 8, "n_r": 6},
        ]

        rendered = markdown_table(rows, ["drug", "n_s", "n_r"])

        self.assertIn("| drug | n_s | n_r |", rendered)
        self.assertIn("| Ciprofloxacin | 10 | 5 |", rendered)
        self.assertIn("| Cefotaxime | 8 | 6 |", rendered)


if __name__ == "__main__":
    unittest.main()
