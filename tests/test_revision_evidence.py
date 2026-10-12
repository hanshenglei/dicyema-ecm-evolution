import importlib.util
import re
import tempfile
import unittest
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('revision_analysis',ROOT/'scripts/publication_revision_analysis.py')
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class RevisionEvidenceTests(unittest.TestCase):
    def test_inclusive_local_motif_intervals(self):
        self.assertEqual(MODULE.intervals('MEME-1:5-12;FG-GAP:20-31'),[('MEME-1',5,12),('FG-GAP',20,31)])

    def test_busco_duplicates_are_one_marker_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'table.tsv'
            path.write_text('# heading\na\tDuplicated\tp1\na\tDuplicated\tp2\nb\tMissing\n')
            self.assertEqual(MODULE.marker_status(path),{'a':'Duplicated','b':'Missing'})

    def test_conflicting_marker_status_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'table.tsv'
            path.write_text('a\tComplete\np\tMissing\na\tFragmented\n')
            with self.assertRaises(ValueError):MODULE.marker_status(path)

    def test_predefined_patterns_are_not_meme_ids(self):
        self.assertIsNotNone(re.fullmatch('HE..H..G..HE','HELMHAIGFYHE'))
        self.assertIsNone(re.fullmatch('HE..H..G..HE','HELMHAIGFYHQ'))
        self.assertIsNotNone(re.fullmatch('NP.Y','NPIY'))


if __name__=='__main__':unittest.main()
