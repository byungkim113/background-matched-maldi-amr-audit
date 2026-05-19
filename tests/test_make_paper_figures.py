import pathlib
import tempfile
import unittest
from unittest import mock

from scripts import make_paper_figures


class MakePaperFiguresTests(unittest.TestCase):
    def test_output_dir_is_forwarded_to_builder(self):
        with tempfile.TemporaryDirectory() as tmp:
            builder = pathlib.Path(tmp) / "builder.py"
            out_dir = pathlib.Path(tmp) / "paper_outputs"

            with mock.patch.object(make_paper_figures.subprocess, "run") as run:
                make_paper_figures.main(
                    [
                        "--builder",
                        str(builder),
                        "--output-dir",
                        str(out_dir),
                    ]
                )

        command = run.call_args.args[0]
        self.assertIn("--output-dir", command)
        self.assertIn(str(out_dir), command)


if __name__ == "__main__":
    unittest.main()
