import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from publication_coordinates import intersection, merged
from publication_analyze import classification
from publication_domain_validation import GROUPS, CORE, MODEL_GROUP


class EvidenceTests(unittest.TestCase):
    def test_interval_union_is_inclusive(self):
        self.assertEqual(merged([(5, 8), (1, 4), (20, 25)]), [(1, 8), (20, 25)])

    def test_intersection_does_not_double_count(self):
        self.assertEqual(intersection([(1, 10), (5, 15)], [(8, 20)]), 8)
        self.assertEqual(intersection([(1, 5)], [(6, 9)]), 0)

    def test_domain_reclassification_uses_evidence(self):
        self.assertEqual(classification(['GCN5L1'])[1], 'intracellular_or_nonECM')
        self.assertEqual(classification(['Integrin_beta'])[1], 'receptor_domain')
        self.assertEqual(classification([])[1], 'unresolved')

    def test_fixed_profile_contract(self):
        self.assertEqual(len(MODEL_GROUP), 40)
        self.assertEqual(len(CORE), 14)
        self.assertTrue(CORE <= set(MODEL_GROUP))


if __name__ == '__main__':
    unittest.main()
