# -*- coding: utf-8 -*-

import logging

from pybel.constants import *
from .manager import Manager

log = logging.getLogger(__name__)


def enrich_rnas(graph, manager=None):
    """Adds all of the miRNA inhibitors of the proteins in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param manager: A mirTarBase database manager
    :type manager: None or str or Manager
    """
    manager = Manager.ensure(manager)

    for node, data in graph.nodes(data=True):
        if data[FUNCTION] != RNA:
            continue

        if NAMESPACE not in data:
            continue

        if data[NAMESPACE] == 'HGNC':
            if IDENTIFIER in data:
                target = manager.query_target_by_hgnc_identifier(data[IDENTIFIER])
            elif NAME in data:
                target = manager.query_target_by_hgnc_symbol(data[NAME])
            else:
                raise IndexError

        elif data[NAMESPACE] in {'EGID', 'ENTREZ'}:
            if IDENTIFIER in data:
                target = manager.query_target_by_entrez_id(data[IDENTIFIER])
            elif NAME in data:
                target = manager.query_target_by_entrez_id(data[NAME])
            else:
                raise IndexError
        else:
            log.warning("Unable to map namespace: %s", data[NAMESPACE])
            continue

        if target is None:
            log.warning("Unable to find RNA: %s:%s", data[NAMESPACE], data[NAME])
            continue

        for interaction in target.interactions:
            mirna_data = interaction.mirna.serialize_to_bel()
            mirna_tuple = graph.add_node_from_data(mirna_data)
            graph.add_edge(mirna_tuple, node, attr_dict={
                RELATION: DIRECTLY_DECREASES,
                EVIDENCE: '...',  # FIXME
                CITATION: {
                    CITATION_TYPE: CITATION_TYPE_PUBMED,
                    CITATION_REFERENCE: interaction.evidence.reference,
                },
                ANNOTATIONS: {
                    'Experiment': interaction.evidence.experiment,
                    'SupportType': interaction.evidence.support,
                }
            })
