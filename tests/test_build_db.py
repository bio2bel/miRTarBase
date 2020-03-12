# -*- coding: utf-8 -*-

"""Tests for Bio2BEL miRTarBase."""

from bio2bel_mirtarbase.models import Evidence, HGNC, MIRBASE_MATURE, Mirna, NCBIGENE, Species, Target
from pybel import BELGraph
from pybel.constants import FUNCTION, IDENTIFIER, NAME, NAMESPACE
from pybel.dsl import BaseAbundance, MicroRna, Rna
from tests.constants import TemporaryFilledCacheMixin

hif1a_symbol = 'HIF1A'

hif1a_hgnc_name = Rna(name=hif1a_symbol, namespace=HGNC)
hif1a_hgnc_identifier = Rna(identifier='4910', namespace=HGNC)
hif1a_entrez_name = Rna(name='3091', namespace=NCBIGENE)
hif1a_entrez_identifier = Rna(identifier='3091', namespace=NCBIGENE)

mi2_data = MicroRna(name='hsa-miR-20a-5p', namespace=MIRBASE_MATURE)
mi5_data = MicroRna(name='mmu-miR-124-3p', namespace=MIRBASE_MATURE)


class TestBuildDatabase(TemporaryFilledCacheMixin):
    """Test the database."""

    def test_count_mirnas(self):
        """Test the number of miRNAs."""
        self.assertEqual(5, self.manager.count_mirnas())

    def test_count_targets(self):
        """Test the number of targets."""
        self.assertEqual(6, self.manager.count_targets())

    def test_count_interactions(self):
        """Test the number of interactions."""
        self.assertEqual(6, self.manager.count_interactions())

    def test_count_evidences(self):
        """Test the number of evidences."""
        self.assertEqual(10, self.manager.count_evidences())

    def test_count_species(self):
        """Test the number of species."""
        self.assertEqual(3, self.manager.session.query(Species).count())

    def test_evidence(self):
        """Test the populate function of the database manager."""
        ev2 = self.manager.session.query(Evidence).filter(Evidence.reference == '18619591').first()
        self.assertIsNotNone(ev2)
        self.assertEqual("Luciferase reporter assay//qRT-PCR//Western blot//Reporter assay;Microarray", ev2.experiment)

    def check_mir5(self, model: Mirna):
        """Help check the model has the right information for mmu-miR-124-3p."""
        self.assertIsNotNone(model)
        self.assertEqual("mmu-miR-124-3p", model.name)
        self.assertTrue(any('MIRT000005' == interaction.mirtarbase_id for interaction in model.interactions))

        bel_data = model.as_bel()

        self.assertEqual(mi5_data.function, bel_data.function)
        self.assertEqual(mi5_data.name, bel_data.name)
        self.assertEqual(mi5_data.namespace, bel_data.namespace)

    def test_mirna_by_mirtarbase_id(self):
        """Test getting an miRNA by a miRTarBase relationship identifier."""
        mi5 = self.manager.query_mirna_by_mirtarbase_identifier('MIRT000005')
        self.check_mir5(mi5)

    def check_mir2(self, model: Mirna):
        """Help check the model has the right information for mmu-miR-124-3p."""
        self.assertIsNotNone(model)
        self.assertEqual("hsa-miR-20a-5p", model.name)
        self.assertEqual(2, len(model.interactions))
        self.assertTrue(any('MIRT000002' == interaction.mirtarbase_id for interaction in model.interactions))

        bel_data = model.as_bel()

        self.assertEqual(mi2_data[FUNCTION], bel_data[FUNCTION])
        self.assertEqual(mi2_data[NAME], bel_data[NAME])
        self.assertEqual(mi2_data[NAMESPACE], bel_data[NAMESPACE])

    def test_mirna_2_by_mirtarbase_id(self):
        """Test getting an miRNA by a miRTarBase relationship identifier."""
        mi2 = self.manager.query_mirna_by_mirtarbase_identifier('MIRT000002')
        self.check_mir2(mi2)

    def test_target(self):
        """Test getting a target by Entrez Gene identifier."""
        target = self.manager.query_target_by_entrez_id('7852')
        self.assertIsNotNone(target)
        self.assertEqual("CXCR4", target.name)
        self.assertIsNotNone(target.hgnc_id)
        self.assertEqual("2561", target.hgnc_id)

    def check_hif1a(self, model: Target):
        """Help check the model has all the right information for HIF1A.

        :type model: Target
        """
        self.assertIsNotNone(model)
        self.assertEqual('HIF1A', model.name)
        self.assertIsNotNone(model.hgnc_id)
        self.assertEqual('4910', model.hgnc_id)
        self.assertIsNotNone(model.entrez_id)
        self.assertEqual('3091', model.entrez_id)

        self.assertEqual(1, len(model.interactions))  # all different evidences to hsa-miR-20a-5p

    def test_target_by_entrez(self):
        """Test getting a target by Entrez Gene identifier."""
        model = self.manager.query_target_by_entrez_id('3091')
        self.check_hif1a(model)

    def test_target_by_hgnc_id(self):
        """Test getting a target by Entrez Gene identifier."""
        model = self.manager.query_target_by_hgnc_identifier('4910')
        self.check_hif1a(model)

    def test_target_by_hgnc_symbol(self):
        """Test getting a target by HGNC symbol."""
        model = self.manager.query_target_by_hgnc_symbol(hif1a_symbol)
        self.check_hif1a(model)

    def help_enrich_hif1a(self, node: BaseAbundance):
        """Help check that different versions of HIF1A can be enriched properly.

        :param pybel.dsl.BaseAbundance node: A PyBEL data dictionary
        """
        self.assertIsInstance(node, BaseAbundance)
        self.assertTrue(NAME in node or IDENTIFIER in node,
                        msg='Node missing information: {}'.format(node))

        graph = BELGraph()
        graph.add_node_from_data(node)
        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

        self.manager.enrich_rnas(graph)  # should enrich with the HIF1A - hsa-miR-20a-5p interaction
        self.assertEqual(2, graph.number_of_nodes(), msg=f"""
        Nodes:
        {", ".join(map(str, graph))}
        """)
        self.assertEqual(3, graph.number_of_edges())

        self.assertIn(mi2_data, graph)
        self.assertTrue(graph.has_edge(mi2_data, node))

    def test_enrich_hgnc_symbol(self):
        """Test enrichment of an HGNC gene symbol node."""
        self.help_enrich_hif1a(hif1a_hgnc_name)

    def test_enrich_hgnc_identifier(self):
        """Test enrichment of an HGNC identifier node."""
        self.help_enrich_hif1a(hif1a_hgnc_identifier)

    def test_enrich_entrez_name(self):
        """Test enrichment of an Entrez Gene node."""
        self.help_enrich_hif1a(hif1a_entrez_name)

    def test_enrich_entrez_id(self):
        """Test enrichment of an Entrez Gene node."""
        self.help_enrich_hif1a(hif1a_entrez_identifier)
