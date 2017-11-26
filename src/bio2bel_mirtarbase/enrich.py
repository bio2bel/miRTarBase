# -*- coding: utf-8 -*-

import logging

from .manager import Manager

log = logging.getLogger(__name__)


def enrich_mirnas(graph, manager=None):
    """Adds all target RNAs to the miRNA nodes in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param manager: A miRTarBase database manager
    :type manager: None or str or Manager
    """
    manager = Manager.ensure(manager)
    manager.enrich_mirnas(graph)


def enrich_rnas(graph, manager=None):
    """Adds all of the miRNA inhibitors of the RNA nodes in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param manager: A miRTarBase database manager
    :type manager: None or str or Manager
    """
    manager = Manager.ensure(manager)
    manager.enrich_rnas(graph)
