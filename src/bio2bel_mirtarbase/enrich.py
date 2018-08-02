# -*- coding: utf-8 -*-

"""Convenience functions for enriching BEL graphs."""

from .manager import Manager

__all__ = [
    'enrich_mirnas',
    'enrich_rnas',
]


def enrich_mirnas(graph, manager=None):
    """Add all target RNAs to the miRNA nodes in the graph.

    :param pybel.BELGraph graph: A BEL graph
    :param manager: A miRTarBase database manager
    :type manager: Optional[Manager]
    """
    if manager is None:
        manager = Manager()
    manager.enrich_mirnas(graph)


def enrich_rnas(graph, manager=None):
    """Add all of the miRNA inhibitors of the RNA nodes in the graph.

    :param pybel.BELGraph graph: A BEL graph
    :param manager: A miRTarBase database manager
    :type manager: Optional[Manager]
    """
    if manager is None:
        manager = Manager()
    manager.enrich_rnas(graph)
