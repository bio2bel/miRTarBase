# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from pyhgnc.manager.database import DbManager
from pyhgnc.manager.query import QueryManager

from bio2bel_mirtarbase.manager import Manager

dir_path = os.path.dirname(os.path.realpath(__file__))
test_xls_path = os.path.join(dir_path, 'test_mirtarbase.xlsx')
test_hgnc_path = os.path.join(dir_path, 'test_hgnc.json')


class HgncManager(DbManager, QueryManager):
    pass


class TemporaryCacheClassMixin(unittest.TestCase):
    """
    :cvar Manager manager: The miRTarBase database manager
    """

    @classmethod
    def setUpClass(cls):
        """Create temporary file"""
        cls.fd, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path
        cls.pyhgnc_manager = HgncManager(connection=cls.connection)
        cls.manager = Manager(connection=cls.connection)

    @classmethod
    def tearDownClass(cls):
        """Closes the connection in the manager and deletes the temporary database"""
        cls.manager.session.close()
        os.close(cls.fd)
        os.remove(cls.path)
