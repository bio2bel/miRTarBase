# -*- coding: utf-8 -*-

"""Manager for Bio2BEL miRTarBase."""

import itertools as itt
import logging
import time
from typing import List, Mapping, Optional

import pandas as pd
from tqdm import tqdm

import pybel.dsl
from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from pybel import BELGraph
from .constants import MIRTARBASE_VERSION, MODULE_NAME
from .models import Base, Evidence, Interaction, Mirna, Species, Target
from .parser import get_preprocessed_df

__all__ = [
    'Manager',
    'get_bel',
]

logger = logging.getLogger(__name__)

VALID_ENTREZ_NAMESPACES = {'egid', 'eg', 'entrez', 'ncbigene'}


def _reverse_dict(d):
    return {v: k for k, v in d.items()}


def get_bel() -> BELGraph:
    """Serialize miRTarBase as BEL."""
    graph = BELGraph(
        name='miRTarBase',
        version=MIRTARBASE_VERSION,
    )
    df = get_preprocessed_df()

    it = df[[
        'mirbase_id', 'mirbase_name', 'mirbase_mature_id', 'mirna_name',
        'target_entrez_id', 'target_name', 'pubmed_id', 'support_type', 'experiments',
    ]]

    for (
        mirbase_id, mirbase_name, mirbase_mature_id, mirbase_mature_name,
        target_entrez_id, target_name, pubmed_id, support, experiments
    ) in tqdm(it.values, total=len(df.index), desc='converting miRTarBase'):
        if any(pd.isna(x) for x in [mirbase_id, mirbase_mature_id]):
            continue
        premature_mirna = pybel.dsl.MicroRna(
            namespace='mirbase',
            identifier=mirbase_id,
            name=mirbase_name,
        )
        mature_mirna = pybel.dsl.MicroRna(
            namespace='mirbase.mature',
            identifier=mirbase_mature_id,
            name=mirbase_mature_name,
        )
        target = pybel.dsl.Rna(
            namespace='ncbigene',
            identifier=target_entrez_id,
            name=target_name,
        )
        graph.add_is_a(mature_mirna, premature_mirna)
        graph.add_directly_decreases(
            mature_mirna, target,
            citation=pubmed_id,
            evidence='From miRTarBase',
            annotations={
                'experiments': set(experiments.split('//')),
                'support': support,
            },
        )

    return graph


class Manager(AbstractManager, BELManagerMixin, FlaskMixin):
    """miRNA-target interactions."""

    module_name = MODULE_NAME
    _base = Base
    edge_model = Interaction
    flask_admin_models = [Mirna, Target, Species, Interaction, Evidence]

    def is_populated(self) -> bool:
        """Check if the database is already populated."""
        return 0 < self.count_mirnas()

    def populate(self, source: Optional[str] = None) -> None:
        """Populate database with the data from miRTarBase.

        :param source: path or link to data source needed for :func:`get_data`
        """
        t = time.time()
        logger.info('getting data')
        df = get_preprocessed_df()
        logger.info('got data in %.2f seconds', time.time() - t)

        # Make species
        taxonomy_id_to_species = {}
        it = set(map(tuple, itt.chain(
            df[['mirna_species_taxonomy_id', 'mirna_species_name']].drop_duplicates().values,
            df[['target_species_taxonomy_id', 'target_species_name']].drop_duplicates().values
        )))
        it = tqdm(it, desc='serializing species')
        for taxonomy_id, name in it:
            species = taxonomy_id_to_species[taxonomy_id] = Species(name=name, taxonomy_id=taxonomy_id)
            self.session.add(species)

        # Make mirna mapping
        mirbase_mature_id_to_mirna = {}
        it = df[[
            'mirbase_mature_id', 'mirna_name', 'mirbase_id', 'mirbase_name', 'mirna_species_taxonomy_id',
        ]].drop_duplicates()
        it = tqdm(it.values, total=len(it.index), desc='serializing micro-rnas')
        for mirbase_mature_id, mirbase_mature_name, mirbase_id, mirbase_name, taxonomy_id in it:
            mirna = mirbase_mature_id_to_mirna[mirbase_mature_id] = Mirna(
                mirbase_mature_name=mirbase_mature_name,
                mirbase_mature_id=mirbase_mature_id,
                mirbase_name=mirbase_name,
                mirbase_id=mirbase_id,
                species=taxonomy_id_to_species[taxonomy_id],
            )
            self.session.add(mirna)

        # make target mapping
        entrez_id_to_target = {}
        it = df[['target_entrez_id', 'target_name', 'target_species_taxonomy_id']].drop_duplicates()
        it = tqdm(it.values, total=len(it.index), desc='serializing targets')
        for entrez_id, name, taxonomy_id in it:
            target = entrez_id_to_target[entrez_id] = Target(
                entrez_id=entrez_id,
                name=name,
                species=taxonomy_id_to_species[taxonomy_id],
            )
            self.session.add(target)

        mirtarbase_id_to_interaction = {}
        it = df[['mirtarbase_id', 'mirbase_mature_id', 'target_entrez_id']].drop_duplicates()
        it = tqdm(it.values, total=len(it.index), desc='serializing interactions')
        for mirtarbase_id, mirbase_mature_id, entrez_id in it:
            interaction = mirtarbase_id_to_interaction[mirtarbase_id] = Interaction(
                mirtarbase_id=mirtarbase_id,
                mirna=mirbase_mature_id_to_mirna[mirbase_mature_id],
                target=entrez_id_to_target[entrez_id],
            )
            self.session.add(interaction)

        t = time.time()
        it = df[['mirtarbase_id', 'experiments', 'support_type', 'pubmed_id']]
        it = tqdm(it.values, total=len(it.index), desc='serializing evidences')
        for mirtarbase_id, experiment, support, pubmed_id in it:
            evidence = Evidence(
                experiment=experiment,
                support=support,
                reference=pubmed_id,
                interaction=mirtarbase_id_to_interaction[mirtarbase_id],
            )
            self.session.add(evidence)

        logger.info('built models in %.2f seconds', time.time() - t)

        logger.info('committing models')
        t = time.time()
        self.session.commit()
        logger.info('committed after %.2f seconds', time.time() - t)

    def count_targets(self) -> int:
        """Count the number of targets in the database."""
        return self._count_model(Target)

    def count_mirnas(self) -> int:
        """Count the number of miRNAs in the database."""
        return self._count_model(Mirna)

    def count_interactions(self) -> int:
        """Count the number of interactions in the database."""
        return self._count_model(Interaction)

    def count_evidences(self) -> int:
        """Count the number of evidences in the database."""
        return self._count_model(Evidence)

    def list_evidences(self) -> List[Evidence]:
        """List the evidences in the database."""
        return self._list_model(Evidence)

    def count_species(self) -> int:
        """Count the number of species in the database."""
        return self._count_model(Species)

    def summarize(self) -> Mapping[str, int]:
        """Return a summary dictionary over the content of the database."""
        return dict(
            targets=self.count_targets(),
            mirnas=self.count_mirnas(),
            species=self.count_species(),
            interactions=self.count_interactions(),
            evidences=self.count_evidences(),
        )

    def query_mirna_by_mirtarbase_identifier(self, mirtarbase_id: str) -> Optional[Mirna]:
        """Get an miRNA by the miRTarBase interaction identifier.

        :param mirtarbase_id: An miRTarBase interaction identifier
        """
        interaction = self.session.query(Interaction).filter(Interaction.mirtarbase_id == mirtarbase_id).one_or_none()
        if interaction is not None:
            return interaction.mirna

    def query_mirna_by_mirtarbase_name(self, name: str) -> Optional[Mirna]:
        """Get an miRNA by its miRTarBase name.

        :param name: An miRTarBase name
        """
        return self.session.query(Mirna).filter(Mirna.name == name).one_or_none()

    def query_mirna_by_hgnc_identifier(self, hgnc_id: str) -> Optional[Mirna]:
        """Query for a miRNA by its HGNC identifier.

        :param hgnc_id: HGNC gene identifier
        """
        raise NotImplementedError

    def query_mirna_by_hgnc_symbol(self, hgnc_symbol: str) -> Optional[Mirna]:
        """Query for a miRNA by its HGNC gene symbol.

        :param hgnc_symbol: HGNC gene symbol
        """
        raise NotImplementedError

    def query_target_by_entrez_id(self, entrez_id: str) -> Optional[Target]:
        """Query for one target.

        :param entrez_id: Entrez gene identifier
        """
        return self.session.query(Target).filter(Target.entrez_id == entrez_id).one_or_none()

    def query_target_by_hgnc_symbol(self, hgnc_symbol: str) -> Optional[Target]:
        """Query for one target.

        :param hgnc_symbol: HGNC gene symbol
        """
        return self.session.query(Target).filter(Target.hgnc_symbol == hgnc_symbol).one_or_none()

    def query_target_by_hgnc_identifier(self, hgnc_id: str) -> Optional[Target]:
        """Query for one target.

        :param hgnc_id: HGNC gene identifier
        """
        return self.session.query(Target).filter(Target.hgnc_id == hgnc_id).one_or_none()

    def _enrich_rna_handle_hgnc(self, identifier, name):
        if identifier:
            return self.query_target_by_hgnc_identifier(identifier)
        if name:
            return self.query_target_by_hgnc_symbol(name)
        raise IndexError

    def _enrich_rna_handle_entrez(self, identifier, name):
        if identifier:
            return self.query_target_by_entrez_id(identifier)
        if name:
            return self.query_target_by_entrez_id(name)
        raise IndexError

    def enrich_rnas(self, graph: BELGraph):
        """Add all of the miRNA inhibitors of the RNA nodes in the graph."""
        logger.debug('enriching inhibitors of RNA')
        count = 0

        for node in list(graph):
            if not isinstance(node, pybel.dsl.Rna):
                continue
            namespace = node.namespace
            if namespace is None:
                continue

            identifier = node.identifier
            name = node.name

            if namespace.lower() == 'hgnc':
                target = self._enrich_rna_handle_hgnc(identifier, name)
            elif namespace.lower() in VALID_ENTREZ_NAMESPACES:
                target = self._enrich_rna_handle_entrez(identifier, name)
            else:
                logger.warning("Unable to map namespace: %s", namespace)
                continue

            if target is None:
                logger.warning("Unable to find RNA: %s:%s", namespace, node)
                continue

            for interaction in target.interactions:
                for evidence in interaction.evidences:
                    count += 1
                    evidence._add_to_graph(graph, evidence.interaction.mirna.as_bel(), node)

        logger.debug('added %d MTIs', count)

    def enrich_mirnas(self, graph: BELGraph):
        """Add all target RNAs to the miRNA nodes in the graph."""
        logger.debug('enriching miRNA targets')
        count = 0

        mirtarbase_names = set()

        for node in graph:
            if not isinstance(node, pybel.dsl.MicroRna):
                continue

            namespace = node.namespace
            name = node.name
            # identifier = node.identifier

            if namespace.lower() == 'mirtarbase':
                if name is not None:
                    mirtarbase_names.add(name)
                raise IndexError('no usable identifier for {}'.format(node))

            elif namespace.lower() in {'mirbase', 'hgnc'} | VALID_ENTREZ_NAMESPACES:
                logger.debug('not yet able to map %s', namespace)
                continue

            else:
                logger.debug("unable to map namespace: %s", namespace)
                continue

        if not mirtarbase_names:
            logger.debug('no mirnas found')
            return

        query = self.get_mirna_interaction_evidences().filter(Mirna.filter_name_in(mirtarbase_names))
        for mirna, interaction, evidence in query:
            count += 1
            evidence.add_to_graph(graph)

        logger.debug('added %d MTIs', count)

    def get_mirna_interaction_evidences(self):
        """Get interaction evidences."""
        return self.session.query(Evidence).join(Evidence.interaction).join(Interaction.mirna)

    def to_bel(self) -> BELGraph:
        """Serialize miRNA-target interactions to BEL."""
        return get_bel()
