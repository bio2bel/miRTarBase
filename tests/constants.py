# -*- coding: utf-8 -*-

"""Test constants for Bio2BEL miRTarBase."""

import json
import os

import bio2bel_hgnc
import bio2bel_mirbase
import bio2bel_mirtarbase
from bio2bel.testing import AbstractTemporaryCacheClassMixin

HERE = os.path.dirname(os.path.realpath(__file__))
TEST_MIRTARBASE_EXCEL = os.path.join(HERE, 'test_mirtarbase.xlsx')
TEST_HGNC_JSON = os.path.join(HERE, 'test_hgnc.json')
TEST_MIRBASE_JSON = os.path.join(HERE, 'test_mirbase.json')


class TemporaryFilledCacheMixin(AbstractTemporaryCacheClassMixin):
    """A test case that holds a temporary database."""

    Manager = bio2bel_mirtarbase.Manager
    manager: Manager
    hgnc_manager: bio2bel_hgnc.Manager
    mirbase_manager: bio2bel_mirbase.Manager

    @classmethod
    def setUpClass(cls):
        """Create temporary file and populate database."""
        super().setUpClass()

        cls.hgnc_manager = bio2bel_hgnc.Manager(connection=cls.connection)
        cls.hgnc_manager._create_tables()
        json_data = cls.hgnc_manager.load_hgnc_json(hgnc_file_path=TEST_HGNC_JSON)
        cls.hgnc_manager.insert_hgnc(hgnc_dict=json_data, silent=True, low_memory=False)

        cls.mirbase_manager = bio2bel_mirbase.Manager(connection=cls.connection)
        with open(TEST_MIRBASE_JSON) as file:
            mirbase_list = json.load(file)
            cls.mirbase_manager._populate_list(mirbase_list)

    @classmethod
    def populate(cls):
        """Fill the HGNC and mirTarBase databases.

        Contents of the test Excel sheet:

        miRTarBase ID	miRNA	Species (miRNA)	Target Gene	Target Gene (Entrez Gene ID)	Species (Target Gene)	Experiments	Support Type	References (PMID)
        MIRT000002	hsa-miR-20a-5p	Homo sapiens	HIF1A	3091	Homo sapiens	Luciferase reporter assay//Western blot//Northern blot//qRT-PCR	Functional MTI	18632605
        MIRT000002	hsa-miR-20a-5p	Homo sapiens	HIF1A	3091	Homo sapiens	Luciferase reporter assay//qRT-PCR//Western blot	Functional MTI	23911400
        MIRT000002	hsa-miR-20a-5p	Homo sapiens	HIF1A	3091	Homo sapiens	HITS-CLIP	Functional MTI (Weak)	22473208
        MIRT000178	hsa-miR-20a-5p	Homo sapiens	TCEAL1	9338	Homo sapiens	Luciferase reporter assay//Microarray//Northern blot//qRT-PCR//Western blot	Functional MTI	23059786
        MIRT000004	dme-miR-8-3p	Drosophila melanogaster	ush	33225	Drosophila melanogaster	qRT-PCR//Luciferase reporter assay//Western blot	Functional MTI	20005803
        MIRT000005	mmu-miR-124-3p	Mus musculus	Itgb1	16412	Mus musculus	Luciferase reporter assay//Microarray//qRT-PCR	Functional MTI	18042700
        MIRT000005	mmu-miR-124-3p	Mus musculus	Itgb1	16412	Mus musculus	Luciferase reporter assay//qRT-PCR//Western blot//Reporter assay;Microarray	Functional MTI	18619591
        MIRT000006	hsa-miR-146a-5p	Homo sapiens	CXCR4	7852	Homo sapiens	qRT-PCR//Luciferase reporter assay//Western blot	Functional MTI	18568019
        MIRT000006	hsa-miR-146a-5p	Homo sapiens	CXCR4	7852	Homo sapiens	Microarray	Functional MTI (Weak)	20375304
        MIRT000012	hsa-miR-122-5p	Homo sapiens	CYP7A1	1581	Homo sapiens	qRT-PCR//Luciferase reporter assay	Functional MTI	20351063
        """
        cls.manager.populate(TEST_MIRTARBASE_EXCEL)
