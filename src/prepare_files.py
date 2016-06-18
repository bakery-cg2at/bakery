#!/usr/bin/env python
"""
Copyright (C) 2014-2016 Jakub Krajniak <jkrajniak@gmail.com>

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
import structures

__doc__ = 'Prepare step of bakery'


def _args():
    parser = argparse.ArgumentParser(
        description='Prepares hybrid coordinate and topology files.',
        add_help=True)

    parser.add_argument('--options', help='XML options file', required=True)

    return parser


def main():
    args = _args().parse_args()
    bck_settings = structures.BackmapperSettings2(args.options)

    bck_settings.prepare_hybrid()


if __name__ == '__main__':
    try:
        import ipdb
        with ipdb.launch_ipdb_on_exception():
            main()
    except ImportError:
        main()
