# -*- coding: utf-8 -*-

"""Test constants for Bio2BEL miRTarBase."""

import os

import bio2bel_hgnc
from bio2bel.testing import AbstractTemporaryCacheClassMixin
from bio2bel_mirtarbase.manager import Manager

dir_path = os.path.dirname(os.path.realpath(__file__))
test_xls_path = os.path.join(dir_path, 'test_mirtarbase.xlsx')
test_hgnc_path = os.path.join(dir_path, 'test_hgnc.json')


class TemporaryFilledCacheMixin(AbstractTemporaryCacheClassMixin):
    """A test case that holds a temporary database."""

    Manager = Manager
    manager = None
    hgnc_manager = None

    @classmethod
    def setUpClass(cls):
        """Create temporary file and populate database."""
        super().setUpClass()

        # fill temporary database with test data
        cls.hgnc_manager = bio2bel_hgnc.Manager(connection=cls.connection)
        cls.hgnc_manager._create_tables()
        json_data = cls.hgnc_manager.load_hgnc_json(hgnc_file_path=test_hgnc_path)
        cls.hgnc_manager.insert_hgnc(hgnc_dict=json_data, silent=True)
        cls.manager.populate(test_xls_path)

    @classmethod
    def tearDownClass(cls):
        """Close the connection in the manager and deletes the temporary database."""
        cls.hgnc_manager.session.close()
        cls.manager.session.close()

        super().tearDownClass()
