# -*- coding: utf-8 -*-

import os
import tempfile
import unittest

from bio2bel_mirtarbase.manager import Manager

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
        cls.manager = Manager(connection=cls.connection)

    @classmethod
    def tearDownClass(cls):
        """Closes the connection in the manager and deletes the temporary database"""
        cls.manager.session.close()
        os.close(cls.fd)
        os.remove(cls.path)
