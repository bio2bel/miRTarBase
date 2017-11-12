# -*- coding: utf-8 -*-

import logging

import click

from bio2bel_mirtarbase.manager import Manager


@click.group()
def main():
    """Bio2BEL miRTarBase"""
    logging.basicConfig(level=logging.INFO)


@main.command()
@click.option('-c', '--connection', help="Custom OLS base url")
def populate(connection):
    """Populates the database"""
    m = Manager(connection=connection)
    m.populate()


@main.command()
@click.option('-c', '--connection', help="Custom OLS base url")
def drop(connection):
    """Drops the database"""
    m = Manager(connection=connection)
    m.populate()


@main.command()
def web():
    """Run the web app"""
    from .web import app
    app.run(host='0.0.0.0', port=5000)


if __name__ == '__main__':
    main()
