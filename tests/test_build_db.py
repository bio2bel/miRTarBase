# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from pybel import BELGraph
from pybel.constants import *

from bio2bel_mirtarbase.enrich import enrich_rnas
from bio2bel_mirtarbase.manager import Manager
from bio2bel_mirtarbase.models import Evidence, Mirna, Target

hif1a_symbol = 'HIF1A'

hif1a = {
    FUNCTION: MIRNA,
    NAMESPACE: 'HGNC',
    NAME: hif1a_symbol
}

t1 = MIRNA, 'MTB', 'MIRT000002'

dir_path = os.path.dirname(os.path.realpath(__file__))
test_xls_path = os.path.join(dir_path, 'test.xlsx')


class TemporaryCacheClassMixin(unittest.TestCase):
    """
    :cvar Manager manager: The miRTarBase database manager
    """

    @classmethod
    def setUpClass(cls):
        """Create temporary file"""
        cls.fd, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path

        # create temporary database
        cls.manager = Manager(connection=cls.connection)
        # fill temporary database with test data
        cls.manager.populate(test_xls_path)

        targets = cls.manager.session.query(Target).all()
        log.error('Targets: %s', [t.to_json() for t in targets])

    @classmethod
    def tearDownClass(cls):
        """Closes the connection in the manager and deletes the temporary database"""
        cls.manager.session.close()
        os.close(cls.fd)
        os.remove(cls.path)


class TestBuildDB(TemporaryCacheClassMixin):
    def test_evidence(self):
        """Test the populate function of the database manager"""
        ev2 = self.manager.session.query(Evidence).filter(Evidence.reference == '18619591').first()
        self.assertIsNotNone(ev2)
        self.assertEqual("Luciferase reporter assay//qRT-PCR//Western blot//Reporter assay;Microarray", ev2.experiment)

    def test_mirna(self):
        mi3 = self.manager.session.query(Mirna).filter(Mirna.mirtarbase_id == "MIRT000005").first()
        self.assertIsNotNone(mi3)
        self.assertEqual("mmu-miR-124-3p", mi3.mir_name)

    def test_target(self):
        targ = self.manager.session.query(Target).filter(Target.entrez_id == '7852').first()
        self.assertIsNotNone(targ)
        self.assertEqual("CXCR4", targ.target_gene)

    def test_target_by_entrez(self):
        hif1a_model = self.manager.query_target_by_entrez_id('3091')
        self.assertIsNotNone(hif1a_model)
        self.assertEqual('HIF1A', hif1a_model.target_gene)

    def test_target_by_hgnc(self):
        hif1a_model = self.manager.query_target_by_hgnc_symbol(hif1a_symbol)
        self.assertIsNotNone(hif1a_model)

    def test_enrich(self):
        g = BELGraph()

        hif1a_tuple = g.add_node_from_data(hif1a)

        self.assertEqual(1, g.number_of_nodes())

        enrich_rnas(g, manager=self.manager)
        self.assertEqual(2, g.number_of_nodes())
        self.assertEqual(3, g.number_of_edges())
        self.assertTrue(g.has_edge(t1, hif1a_tuple))


if __name__ == '__main__':
    unittest.main()
