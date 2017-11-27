# -*- coding: utf-8 -*-

import unittest

from bio2bel_mirtarbase.enrich import enrich_rnas
from bio2bel_mirtarbase.manager import build_entrez_map
from bio2bel_mirtarbase.models import Evidence, Interaction, Mirna, Species, Target
from pybel import BELGraph
from pybel.constants import FUNCTION, IDENTIFIER, NAME, NAMESPACE
from pybel.dsl import *
from pybel.parser.canonicalize import node_to_tuple
from tests.constants import TemporaryCacheClassMixin, test_hgnc_path, test_xls_path

hif1a_symbol = 'HIF1A'

hif1a_hgnc_name = rna(name=hif1a_symbol, namespace='HGNC')
hif1a_hgnc_identifier = rna(identifier='4910', namespace='HGNC')
hif1a_entrez_name = rna(name='3091', namespace='EGID')
hif1a_entrez_identifier = rna(identifier='3091', namespace='ENTREZ')
mi2_data = mirna(name='hsa-miR-20a-5p', namespace='MIRTARBASE', identifier='MIRT000002')
mi5_data = mirna(name='mmu-miR-124-3p', namespace='MIRTARBASE', identifier='MIRT000005')


class TemporaryFilledCacheMixin(TemporaryCacheClassMixin):
    """
    :cvar bio2bel_mirtarbase.manager.Manager manager: The miRTarBase database manager
    """

    @classmethod
    def setUpClass(cls):
        """Create temporary file and populate database"""
        super(TemporaryFilledCacheMixin, cls).setUpClass()
        # fill temporary database with test data
        cls.pyhgnc_manager._create_tables()
        json_data = cls.pyhgnc_manager.load_hgnc_json(hgnc_file_path=test_hgnc_path)
        cls.pyhgnc_manager.insert_hgnc(hgnc_dict=json_data, silent=True)

        cls.manager.populate(test_xls_path, pyhgnc_connection=cls.connection)


class TestBuildDatabase(TemporaryFilledCacheMixin):
    def test_count_mirnas(self):
        self.assertEqual(5, self.manager.count_mirnas())

    def test_count_targets(self):
        self.assertEqual(6, self.manager.count_targets())

    def test_count_interactions(self):
        self.assertEqual(6, self.manager.count_interactions())

    def test_count_evidences(self):
        self.assertEqual(10, self.manager.count_evidences())

    def test_count_species(self):
        self.assertEqual(3, self.manager.session.query(Species).count())

    def test_count_hgnc(self):
        self.assertEqual(2, len(self.pyhgnc_manager.hgnc()))

    def get_cxcr4_by_entrez(self):
        model = self.pyhgnc_manager.hgnc(entrez='7852')
        self.assertIsNotNone(model)
        self.assertEqual('CXCR4', model.hgnc_symbol)
        self.assertEqual('7852', model.entrez)

    def get_hif1a_by_entrez(self):
        model = self.pyhgnc_manager.hgnc(entrez='3091')
        self.assertIsNotNone(model)
        self.assertEqual('HIF1A', model.hgnc_symbol)
        self.assertEqual('3091', model.entrez)

    def test_build_map(self):
        emap = build_entrez_map(pyhgnc_connection=self.pyhgnc_manager)
        self.assertEqual(2, len(emap))
        self.assertIn('7852', emap)
        self.assertIn('3091', emap)

    def test_evidence(self):
        """Test the populate function of the database manager"""
        ev2 = self.manager.session.query(Evidence).filter(Evidence.reference == '18619591').first()
        self.assertIsNotNone(ev2)
        self.assertEqual("Luciferase reporter assay//qRT-PCR//Western blot//Reporter assay;Microarray", ev2.experiment)

    def check_mir5(self, model):
        """Checks the model has the right information for mmu-miR-124-3p

        :type model: Mirna
        """
        self.assertIsNotNone(model)
        self.assertEqual("mmu-miR-124-3p", model.name)
        self.assertTrue(any('MIRT000005' == interaction.mirtarbase_id for interaction in model.interactions))

        bel_data = model.serialize_to_bel()

        self.assertEqual(mi5_data[FUNCTION], bel_data[FUNCTION])
        self.assertEqual(mi5_data[NAME], bel_data[NAME])
        self.assertEqual(mi5_data[NAMESPACE], bel_data[NAMESPACE])

    def test_mirna_by_mirtarbase_id(self):
        mi5 = self.manager.query_mirna_by_mirtarbase_identifier('MIRT000005')
        self.check_mir5(mi5)

    def check_mir2(self, model):
        """Checks the model has the right information for mmu-miR-124-3p

        :type model: Mirna
        """
        self.assertIsNotNone(model)
        self.assertEqual("hsa-miR-20a-5p", model.name)
        self.assertEqual(2, len(model.interactions))
        self.assertTrue(any('MIRT000002' == interaction.mirtarbase_id for interaction in model.interactions))

        bel_data = model.serialize_to_bel()

        self.assertEqual(mi2_data[FUNCTION], bel_data[FUNCTION])
        self.assertEqual(mi2_data[NAME], bel_data[NAME])
        self.assertEqual(mi2_data[NAMESPACE], bel_data[NAMESPACE])

    def test_mirna_2_by_mirtarbase_id(self):
        mi2 = self.manager.query_mirna_by_mirtarbase_identifier('MIRT000002')
        self.check_mir2(mi2)

    def test_target(self):
        target = self.manager.query_target_by_entrez_id('7852')
        self.assertIsNotNone(target)
        self.assertEqual("CXCR4", target.name)
        self.assertEqual("2561", target.hgnc_id)

    def check_hif1a(self, model):
        """Checks the model has all the right information for HIF1A

        :type model: Target
        """
        self.assertIsNotNone(model)
        self.assertEqual('HIF1A', model.name)
        self.assertEqual('4910', model.hgnc_id)
        self.assertEqual('HIF1A', model.hgnc_symbol)
        self.assertEqual('3091', model.entrez_id)

        self.assertEqual(1, len(model.interactions))  # all different evidences to hsa-miR-20a-5p

    def test_target_by_entrez(self):
        model = self.manager.query_target_by_entrez_id('3091')
        self.check_hif1a(model)

    def test_target_by_hgnc_id(self):
        model = self.manager.query_target_by_hgnc_identifier('4910')
        self.check_hif1a(model)

    def test_target_by_hgnc_symbol(self):
        model = self.manager.query_target_by_hgnc_symbol(hif1a_symbol)
        self.check_hif1a(model)

    def help_enrich_hif1a(self, node_data):
        """Checks that different versions of HIF1A can be enriched properly

        :param dict node_data: A PyBEL data dictionary
        """
        self.assertTrue(NAME in node_data or IDENTIFIER in node_data,
                        msg='Node missing information: {}'.format(node_data))

        graph = BELGraph()

        hif1a_tuple = graph.add_node_from_data(node_data)
        self.assertEqual(1, graph.number_of_nodes())

        enrich_rnas(graph, manager=self.manager)  # should enrich with the HIF1A - hsa-miR-20a-5p interaction
        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(3, graph.number_of_edges())

        self.assertTrue(graph.has_node_with_data(mi2_data))

        t1_tuple = node_to_tuple(mi2_data)
        self.assertTrue(graph.has_edge(t1_tuple, hif1a_tuple))

    def test_enrich_hgnc_symbol(self):
        self.help_enrich_hif1a(hif1a_hgnc_name)

    def test_enrich_hgnc_identifier(self):
        self.help_enrich_hif1a(hif1a_hgnc_identifier)

    def test_enrich_entrez_name(self):
        self.help_enrich_hif1a(hif1a_entrez_name)

    def test_enrich_entrez_id(self):
        self.help_enrich_hif1a(hif1a_entrez_identifier)


if __name__ == '__main__':
    unittest.main()
