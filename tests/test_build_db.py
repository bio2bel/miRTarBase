# -*- coding: utf-8 -*-

import unittest

from bio2bel_mirtarbase.enrich import enrich_rnas
from bio2bel_mirtarbase.models import Evidence, Interaction, Mirna, Species, Target
from pybel import BELGraph
from pybel.dsl import *
from pybel.parser.canonicalize import node_to_tuple
from tests.constants import TemporaryCacheClassMixin, test_xls_path

hif1a_symbol = 'HIF1A'

hif1a = rna(name=hif1a_symbol, namespace='HGNC')
t1 = mirna(name='hsa-miR-20a-5p', namespace='MIRTARBASE', identifier='MIRT000002')


class TemporaryFilledCacheMixin(TemporaryCacheClassMixin):
    """
    :cvar Manager manager: The miRTarBase database manager
    """

    @classmethod
    def setUpClass(cls):
        """Create temporary file and populate database"""
        super(TemporaryFilledCacheMixin, cls).setUpClass()
        # fill temporary database with test data
        cls.manager.populate(test_xls_path)


class TestBuildDB(TemporaryFilledCacheMixin):
    def test_count_mirnas(self):
        self.assertEqual(5, self.manager.session.query(Mirna).count())

    def test_count_targets(self):
        self.assertEqual(5, self.manager.session.query(Target).count())

    def test_count_interactions(self):
        self.assertEqual(9, self.manager.session.query(Interaction).count())

    def test_count_species(self):
        self.assertEqual(3, self.manager.session.query(Species).count())

    def test_evidence(self):
        """Test the populate function of the database manager"""
        ev2 = self.manager.session.query(Evidence).filter(Evidence.reference == '18619591').first()
        self.assertIsNotNone(ev2)
        self.assertEqual("Luciferase reporter assay//qRT-PCR//Western blot//Reporter assay;Microarray", ev2.experiment)

    def test_mirna(self):
        mi3 = self.manager.session.query(Mirna).filter(Mirna.mirtarbase_id == "MIRT000005").first()
        self.assertIsNotNone(mi3)
        self.assertEqual("mmu-miR-124-3p", mi3.mirtarbase_name)

    def test_target(self):
        targ = self.manager.session.query(Target).filter(Target.entrez_id == '7852').first()
        self.assertIsNotNone(targ)
        self.assertEqual("CXCR4", targ.target_gene)

    def check_hif1a(self, model):
        """Checks the model has all the right information for HIF1A

        :type model: Target
        """
        self.assertIsNotNone(model)
        self.assertEqual('HIF1A', model.target_gene)
        self.assertEqual('3091', model.entrez_id)

    def test_target_by_entrez(self):
        model = self.manager.query_target_by_entrez_id('3091')
        self.check_hif1a(model)

    def test_target_by_hgnc(self):
        model = self.manager.query_target_by_hgnc_symbol(hif1a_symbol)
        self.check_hif1a(model)

    def test_enrich_hgnc_symbol(self):
        g = BELGraph()

        hif1a_tuple = g.add_node_from_data(hif1a)

        self.assertEqual(1, g.number_of_nodes())

        enrich_rnas(g, manager=self.manager)
        self.assertEqual(2, g.number_of_nodes())
        self.assertEqual(3, g.number_of_edges())

        self.assertTrue(g.has_node_with_data(t1))

        t1_tuple = node_to_tuple(t1)
        self.assertTrue(g.has_edge(t1_tuple, hif1a_tuple))

    def test_enrich_hgnc_id(self):
        pass

    def test_enrich_entrez_id(self):
        pass


if __name__ == '__main__':
    unittest.main()
