# -*- coding: utf-8 -*-

import logging

import click

from bio2bel_mirtarbase.constants import DEFAULT_CACHE_CONNECTION
from bio2bel_mirtarbase.manager import Manager

log = logging.getLogger(__name__)


@click.group()
@click.option('-c', '--connection', help='Defaults to {}'.format(DEFAULT_CACHE_CONNECTION))
@click.pass_context
def main(ctx):
    """miRTarBase to BEL"""
    logging.basicConfig(level=logging.INFO)
    log.setLevel(logging.INFO)
    ctx.obj = Manager(connection=connection)


@main.command()
@click.pass_obj
def populate(manager):
    """Populates the database"""
    manager.populate()


@main.command()
@click.pass_obj
def drop(manager):
    """Drops the database"""
    manager.drop_all()


@main.command()
@click.pass_obj
def summarize(manager):
    """Drops the database"""
    click.echo('miRNAs: {}'.format(manager.count_mirnas()))
    click.echo('Targets: {}'.format(manager.count_targets()))
    click.echo('Species: {}'.format(manager.count_species()))
    click.echo('Interactions: {}'.format(manager.count_interactions()))
    click.echo('Evidences: {}'.format(manager.count_evidences()))

@main.command()
@click.option('-v', '--debug')
@click.option('-p', '--port')
@click.option('-h', '--host')
@click.pass_obj
def web(manager, debug, port, host):
    """Run the admin interface"""
    from .web import get_app
    app = get_app(connection=manager)
    app.run(debug=debug, port=port, host=host)


if __name__ == '__main__':
    main()
