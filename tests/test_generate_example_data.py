import pathlib
import tempfile
import unittest

from scripts import generate_example_data


class GenerateExampleDataTests(unittest.TestCase):
    def test_generated_csv_uses_lf_line_endings(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_out = generate_example_data.OUT
            try:
                generate_example_data.OUT = pathlib.Path(tmp) / "example_predictions.csv"
                generate_example_data.main()
                data = generate_example_data.OUT.read_bytes()
            finally:
                generate_example_data.OUT = old_out

        self.assertNotIn(b"\r\n", data)
        self.assertIn(b"\n", data)


if __name__ == "__main__":
    unittest.main()
