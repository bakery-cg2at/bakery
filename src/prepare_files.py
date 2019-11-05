#!/usr/bin/env python2
"""
Copyright (C) 2014-2017 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of Backmapper.

Backmapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import logging

import structures
from logger import logger

__doc__ = 'Prepare step of bakery'


def _args():
    parser = argparse.ArgumentParser(
        description='Prepares hybrid coordinate and topology files.',
        add_help=True)

    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--options', help='XML options file', required=True)
    parser.add_argument('--allow-no-bonds', help='Allow to skip AT bonds', action='store_true', dest='allow_no_bonds')
    parser.add_argument('--generate-only-graph', help='Generate only NX graph', action='store_true', dest='generate_only_graph')

    return parser


def main():
    args = _args().parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.ERROR)

    bck_settings = structures.BackmapperSettings2(args.options, args.allow_no_bonds, args.generate_only_graph)

    bck_settings.prepare_hybrid()


if __name__ == '__main__':
    main()
