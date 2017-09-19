# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import RNA, MIRNA

from bio2bel_mirtarbase.enrich import enrich_rnas

hif1a_tuple = RNA, 'HGNC', 'HIF1A'
t1 = MIRNA, 'MTB', 'MIRT000002'


class TestEnrich(unittest.TestCase):
    def test_enrich(self):
        g = BELGraph()

        g.add_simple_node(*hif1a_tuple)

        self.assertEqual(1, g.number_of_nodes())

        enrich_rnas(g)
        self.assertEqual(2, g.number_of_nodes())
        self.assertEqual(3, g.number_of_edges())
        self.assertTrue(g.has_edge(t1, hif1a_tuple))
