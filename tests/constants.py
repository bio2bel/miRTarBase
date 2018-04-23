# -*- coding: utf-8 -*-

import os

import bio2bel_hgnc
from bio2bel.testing import make_temporary_cache_class_mixin
from bio2bel_mirtarbase.manager import Manager

dir_path = os.path.dirname(os.path.realpath(__file__))
test_xls_path = os.path.join(dir_path, 'test_mirtarbase.xlsx')
test_hgnc_path = os.path.join(dir_path, 'test_hgnc.json')

_TemporaryCacheClassMixinParent = make_temporary_cache_class_mixin(Manager)


class TemporaryCacheClassMixin(_TemporaryCacheClassMixinParent):
    """
    :cvar Manager manager: The miRTarBase database manager
    """

    @classmethod
    def setUpClass(cls):
        """Create temporary file"""
        super(TemporaryCacheClassMixin, cls).setUpClass()

        cls.hgnc_manager = bio2bel_hgnc.Manager(connection=cls.connection)

    @classmethod
    def tearDownClass(cls):
        """Closes the connection in the manager and deletes the temporary database"""
        cls.hgnc_manager.session.close()
        cls.manager.session.close()

        super(TemporaryCacheClassMixin, cls).tearDownClass()
