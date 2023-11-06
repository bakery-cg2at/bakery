import click
import logging

from src.bakery import structures
from src.bakery.logger import logger


@click.group()
@click.option("--verbose", help="Verbose output", is_flag=True, default=False)
def cli(verbose):
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.ERROR)


@cli.command(help="Prepare configuration files for backmapping")
@click.option("--options", help="XML options file", required=True)
@click.option("--allow-no-bonds", help="Allow to skip AT bonds", is_flag=True, default=False)
@click.option("--generate-only-graph", help="Generate only NetworkX graphs", is_flag=True, default=False)
def prepare_files(options, allow_no_bonds, generate_only_graph):
    bck_settings = structures.BackmapperSettings2(options, allow_no_bonds, generate_only_graph)

    bck_settings.prepare_hybrid()


if __name__ == '__main__':
    cli()
