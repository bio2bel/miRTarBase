# -*- coding: utf-8 -*-

import logging

from pybel.constants import *
from .manager import Manager

log = logging.getLogger(__name__)


def get_name(data):
    if NAME in data:
        return data[NAME]
    elif IDENTIFIER in data:
        return data[IDENTIFIER]


def enrich_rnas(graph, manager=None):
    """Adds all of the miRNA inhibitors of the proteins in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param manager: A miRTarBase database manager
    :type manager: None or str or Manager
    """
    manager = Manager.ensure(manager)

    for node, data in graph.nodes(data=True):
        if data[FUNCTION] != RNA:
            continue

        if NAMESPACE not in data:
            continue

        namespace = data[NAMESPACE]

        if namespace == 'HGNC':
            if IDENTIFIER in data:
                target = manager.query_target_by_hgnc_identifier(data[IDENTIFIER])
            elif NAME in data:
                target = manager.query_target_by_hgnc_symbol(data[NAME])
            else:
                raise IndexError

        elif namespace in {'EGID', 'ENTREZ'}:
            if IDENTIFIER in data:
                target = manager.query_target_by_entrez_id(data[IDENTIFIER])
            elif NAME in data:
                target = manager.query_target_by_entrez_id(data[NAME])
            else:
                raise IndexError
        else:
            log.warning("Unable to map namespace: %s", namespace)
            continue

        if target is None:
            log.warning("Unable to find RNA: %s:%s", namespace, get_name(data))
            continue

        for interaction in target.interactions:
            graph.add_qualified_edge(
                interaction.mirna.serialize_to_bel(),
                node,
                relation=DIRECTLY_DECREASES,
                evidence=interaction.evidence.support,
                citation=str(interaction.evidence.reference),
                annotations={
                    'Experiment': interaction.evidence.experiment,
                    'SupportType': interaction.evidence.support,
                }
            )
