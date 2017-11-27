# -*- coding: utf-8 -*-

import logging

import click

from bio2bel_mirtarbase.constants import DEFAULT_CACHE_CONNECTION
from bio2bel_mirtarbase.manager import Manager

log = logging.getLogger(__name__)


@click.group()
def main():
    """miRTarBase to BEL"""
    logging.basicConfig(level=logging.INFO)
    log.setLevel(logging.INFO)


@main.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
def populate(connection):
    """Populates the database"""
    m = Manager(connection=connection)
    m.populate()


@main.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
def drop(connection):
    """Drops the database"""
    m = Manager(connection=connection)
    m.drop_all()


@main.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
def summarize(connection):
    """Drops the database"""
    m = Manager(connection=connection)
    click.echo('miRNAs: {}'.format(m.count_mirnas()))
    click.echo('Targets: {}'.format(m.count_targets()))
    click.echo('Species: {}'.format(m.count_species()))
    click.echo('Interactions: {}'.format(m.count_interactions()))
    click.echo('Evidences: {}'.format(m.count_evidences()))

@main.command()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--debug')
@click.option('-p', '--port')
@click.option('-h', '--host')
def web(connection, debug, port, host):
    """Run the admin interface"""
    from .web import get_app
    app = get_app(connection=connection)
    app.run(debug=debug, port=port, host=host)


if __name__ == '__main__':
    main()
