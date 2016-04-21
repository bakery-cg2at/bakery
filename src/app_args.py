"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import tools as general_tools


def _args_backmapping():
    parser = general_tools.MyArgParser(
        description='Starts MD simulation and run backmapping',
        fromfile_prefix_chars='@')
    parser.add_argument('--conf', required=True, help='Coordinate file')
    parser.add_argument('--top', '--topology', required=True, help='Topology file',
                        dest='top')
    parser.add_argument('--coord', help='Input h5md coordinate file')
    parser.add_argument('--coord_frame', default=-1, type=int,
                        help='Time frame of the input coordinate h5md file.')
    parser.add_argument('--coord_h5md_group', help='H5MD atoms group',
                        default='atoms')
    parser.add_argument('--node_grid', help='Node grid configuration')
    parser.add_argument('--cell_grid', help='Cell grid configuration')
    parser.add_argument('--skin', type=float, help='Skin value for VerletList')
    parser.add_argument('--initial_resolution',
                        default=0.0, type=float, help='Initial resolution')
    parser.add_argument('--alpha',
                        type=float, default=0.1, help='Alpha parameter')
    parser.add_argument('--eq', type=int, default=10000, help='CG run time')
    parser.add_argument('--thermostat_gamma',
                        type=float, default=0.5,
                        help='Thermostat coupling constant')
    parser.add_argument('--int_step',
                        default=1000, type=int, help='Number of steps in integration phase')
    parser.add_argument('--long',
                        default=10000, type=int, help='Number of steps after backmapping')
    parser.add_argument('--rng_seed',
                        default=-1, type=int, help='Seed for random number generator')
    parser.add_argument('--output_prefix', default='sim', type=str,
                        help='Prefix for output files')
    parser.add_argument('--thermostat',
                        default='lv', choices=('lv', 'vr'),
                        help='Thermostat to use, lv: Langevine, vr: Stochastic velocity rescale')
    parser.add_argument('--temperature',
                        default=423.0, type=float,
                        help='Temperature')
    parser.add_argument('--dt', default=0.001, type=float, help='Integrator time step')
    parser.add_argument('--lj_cutoff',
                        default=1.2, help='Cutoff of atomistic non-bonded interactions',
                        type=float)
    parser.add_argument('--cg_cutoff',
                        default=1.4, help='Cutoff of coarse-grained non-bonded interactions',
                        type=float)
    parser.add_argument('--coulomb_epsilon1', default=1.0, type=float,
                        help='Epsilon 1 parameter for Coulomb reaction field potential')
    parser.add_argument('--coulomb_epsilon2', default=80.0, type=float,
                        help='Epsilon 2 parameter for Coulomb reaction field potential')
    parser.add_argument('--coulomb_kappa', default=1.0, type=float,
                        help='Kappa paramter for coulomb interactions')
    parser.add_argument('--energy_collect', default=1000, type=int,
                        help='Collect energy every (step)')
    parser.add_argument('--energy_collect_bck', default=1, type=int,
                        help='Collect energy every (step) during backmapping')
    parser.add_argument('--trj_collect', default=1000, type=int,
                        help='Collect trajectory every (step)')

    return parser
