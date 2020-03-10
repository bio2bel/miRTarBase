# -*- coding: utf-8 -*-

"""Manager for Bio2BEL miRTarBase."""

import logging
import time
from typing import List, Mapping, Optional

from tqdm import tqdm

import bio2bel_entrez
import bio2bel_hgnc
import bio2bel_mirbase
from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from bio2bel_hgnc.models import HumanGene
from pybel import BELGraph
from pybel.constants import FUNCTION, IDENTIFIER, MIRNA, NAME, NAMESPACE, RNA
from .constants import MIRTARBASE_VERSION, MODULE_NAME
from .models import Base, Evidence, Interaction, Mirna, Species, Target
from .parser import get_data

__all__ = [
    'Manager',
]

logger = logging.getLogger(__name__)

VALID_ENTREZ_NAMESPACES = {'egid', 'eg', 'entrez', 'ncbigene'}


def _build_entrez_map(hgnc_manager: bio2bel_hgnc.Manager) -> Mapping[str, HumanGene]:
    """Build a mapping from entrez gene identifiers to their database models from :py:mod:`bio2bel_hgnc.models`."""
    logger.info('getting entrez mapping')

    t = time.time()
    emap = {
        model.entrez: model
        for model in hgnc_manager.hgnc()
        if model.entrez
    }
    logger.info('got entrez mapping in %.2f seconds', time.time() - t)
    return emap


def _get_name(data):
    if NAME in data:
        return data[NAME]
    elif IDENTIFIER in data:
        return data[IDENTIFIER]


class Manager(AbstractManager, BELManagerMixin, FlaskMixin):
    """miRNA-target interactions."""

    module_name = MODULE_NAME
    _base = Base
    edge_model = Interaction
    flask_admin_models = [Mirna, Target, Species, Interaction, Evidence]

    def is_populated(self) -> bool:
        """Check if the database is already populated."""
        return 0 < self.count_mirnas()

    def populate(self, source: Optional[str] = None, update: bool = False) -> None:
        """Populate database with the data from miRTarBase.

        :param source: path or link to data source needed for :func:`get_data`
        :param update: Should HGNC an miRBase be updated?
        """
        hgnc_manager = bio2bel_hgnc.Manager(connection=self.connection)
        if not hgnc_manager.is_populated() or update:
            hgnc_manager.populate()

        mirbase_manager = bio2bel_mirbase.Manager(connection=self.connection)
        if not mirbase_manager.is_populated() or update:
            mirbase_manager.populate()

        t = time.time()
        logger.info('getting data')
        df = get_data(source)
        logger.info('got data in %.2f seconds', time.time() - t)

        name_mirna = {}
        target_set = {}
        species_set = {}
        interaction_set = {}

        emap = _build_entrez_map(hgnc_manager)

        logger.info('building models')

        t = time.time()
        it = tqdm(df.values, total=len(df.index))
        for mirtarbase_id, mirna_name, mirna_species, gene_name, entrez_id, target_species, exp, sup_type, pubmed in it:
            # create new miRNA instance
            entrez_id = str(int(entrez_id))

            interaction_key = (mirna_name, entrez_id)
            interaction = interaction_set.get(interaction_key)

            if interaction is None:
                mirna = name_mirna.get(mirna_name)
                if mirna is None:
                    species = species_set.get(mirna_species)
                    if species is None:
                        species = species_set[mirna_species] = Species(name=mirna_species)
                        self.session.add(species)

                    mirna = name_mirna[mirna_name] = Mirna(
                        name=mirna_name,
                        species=species,
                    )
                    self.session.add(mirna)

                target = target_set.get(entrez_id)
                if target is None:
                    species = species_set.get(target_species)

                    if species is None:
                        species = species_set[target_species] = Species(name=target_species)
                        self.session.add(species)

                    target = target_set[entrez_id] = Target(
                        entrez_id=entrez_id,
                        species=species,
                        name=gene_name,
                    )

                    if entrez_id in emap:
                        g_first = emap[entrez_id]
                        target.hgnc_symbol = g_first.symbol
                        target.hgnc_id = str(g_first.identifier)

                    self.session.add(target)

                # create new interaction instance
                interaction = interaction_set[interaction_key] = Interaction(
                    mirtarbase_id=mirtarbase_id,
                    mirna=mirna,
                    target=target
                )
                self.session.add(interaction)

            # create new evidence instance
            new_evidence = Evidence(
                experiment=exp,
                support=sup_type,
                reference=pubmed,
                interaction=interaction,
            )
            self.session.add(new_evidence)

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
            if node[FUNCTION] != RNA:
                continue

            namespace = node.get(NAMESPACE)
            if namespace is None:
                continue

            identifier = node.get(IDENTIFIER)
            name = node.get(NAME)

            if namespace.lower() == 'hgnc':
                target = self._enrich_rna_handle_hgnc(identifier, name)
            elif namespace.lower() in VALID_ENTREZ_NAMESPACES:
                target = self._enrich_rna_handle_entrez(identifier, name)
            else:
                logger.warning("Unable to map namespace: %s", namespace)
                continue

            if target is None:
                logger.warning("Unable to find RNA: %s:%s", namespace, _get_name(node))
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
            if node[FUNCTION] != MIRNA or NAMESPACE not in node:
                continue

            namespace = node[NAMESPACE]

            if namespace.lower() == 'mirtarbase':
                if NAME in node:
                    mirtarbase_names.add(node[NAME])
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
        graph = BELGraph(
            name='miRTarBase',
            version=MIRTARBASE_VERSION,
        )
        graph.namespace_pattern.update(dict(
            hgnc=bio2bel_hgnc.Manager.identifiers_pattern,
            entrez=bio2bel_entrez.Manager.identifiers_pattern,
            mirbase=bio2bel_mirbase.Manager.identifiers_pattern,
        ))

        # TODO check if entrez has all species uploaded and optionally populate remaining species
        it = tqdm(
            self.get_mirna_interaction_evidences(),
            total=self.count_evidences(),
            desc='Mapping miRNA-target interactions to BEL',
        )
        for evidence in it:
            evidence.add_to_graph(graph)

        return graph
