miRTarBase to BEL |build| |coverage| |documentation|
====================================================
Serializes text-mined miRNA-target interactions to BEL

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
``bio2bel_mirtarbase`` can be installed easily from `PyPI <https://pypi.python.org/pypi/bio2bel_mirtarbase>`_ with the
following code in your favorite terminal:

.. code-block:: sh

    $ python3 -m pip install bio2bel_mirtarbase

or from the latest code on `GitHub <https://github.com/bio2bel/mirtarbase>`_ with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/bio2bel/mirtarbase.git@master

Citation
--------
Chih-Hung Chou, et al; miRTarBase 2016: updates to the experimentally validated miRNA-target interactions database.
Nucleic Acids Res 2016; 44 (D1): D239-D247. doi: 10.1093/nar/gkv1258

Abstract
--------
As a database, miRTarBase has accumulated more than three hundred and sixty thousand miRNA-target interactions (MTIs),
which are collected by manually surveying pertinent literature after NLP of the text systematically to filter research
articles related to functional studies of miRNAs. Generally, the collected MTIs are validated experimentally by reporter
assay, western blot, microarray and next-generation sequencing experiments. While containing the largest amount of
validated MTIs, the miRTarBase provides the most updated collection by comparing with other similar, previously
developed databases.

Links
-----
- http://mirtarbase.mbc.nctu.edu.tw/


.. |build| image:: https://travis-ci.org/bio2bel/mirtarbase.svg?branch=master
    :target: https://travis-ci.org/bio2bel/mirtarbase
    :alt: Build Status

.. |coverage| image:: https://codecov.io/gh/bio2bel/mirtarbase/coverage.svg?branch=master
    :target: https://codecov.io/gh/bio2bel/mirtarbase?branch=master
    :alt: Coverage Status

.. |documentation| image:: https://readthedocs.org/projects/mirtarbase/badge/?version=latest
    :target: http://mirtarbase.readthedocs.io
    :alt: Documentation Status

.. |climate| image:: https://codeclimate.com/github/bio2bel/mirtarbase/badges/gpa.svg
    :target: https://codeclimate.com/github/bio2bel/mirtarbase
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/bio2bel_mirtarbase.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/bio2bel_mirtarbase.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/bio2bel_mirtarbase.svg
    :alt: MIT License
