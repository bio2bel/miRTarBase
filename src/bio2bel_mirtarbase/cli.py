# -*- coding: utf-8 -*-

import logging

import click


@click.group()
def main():
    """Bio2BEL miRTarBase"""
    logging.basicConfig(level=logging.INFO)


@main.command()
@click.option('-c', '--connection')
def populate(connection):
    """Populates the local miRTarBase database"""
    from bio2bel_mirtarbase.manager import Manager

    manager = Manager(connection=connection)

    manager.populate()
