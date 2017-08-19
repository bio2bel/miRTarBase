# -*- coding: utf-8 -*-

from pybel_tools import pipeline


@pipeline.in_place_mutator
def enrich_proteins(graph):
    """Adds all of the miRNA inhibitors of the proteins in the graph

    :param pybel.BELGraph graph: A BEL graph
    """
    raise NotImplementedError
