# -*- coding: utf-8 -*-

import unittest

from bio2bel_mirtarbase.enrich import enrich_proteins
from pybel import BELGraph
from pybel.constants import PROTEIN, MIRNA

c = PROTEIN, 'HGNC', 'HIF1A'
t1 = MIRNA, 'MTB', 'MIRT000002'


class TestEnrich(unittest.TestCase):
    def test_enrich(self):
        g = BELGraph()

        g.add_simple_node(*c)

        self.assertEqual(1, g.number_of_nodes())

        enrich_proteins(g)

        self.assertTrue(g.has_edge(t1, c))
