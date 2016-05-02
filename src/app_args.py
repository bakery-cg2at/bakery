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
import random
import ast


def _args_md():
    parser = general_tools.MyArgParser(description='Runs classical MD simulation',
                                       fromfile_prefix_chars='@')
    parser.add_argument('--conf', required=True, help='Input .gro coordinate file')
    parser.add_argument('--top', '--topology', required=True, help='Topology file',
                        dest='top')
    parser.add_argument('--node_grid')
    parser.add_argument('--skin', type=float, default=0.16,
                        help='Skin value for Verlet list')
    parser.add_argument('--coord', help='Input coordinate h5md file')
    parser.add_argument('--coord_frame', default=-1, type=int,
                        help='Time frame of input coordinate h5md file')
    parser.add_argument('--run', type=int, default=10000,
                        help='Number of simulation steps')
    parser.add_argument('--int_step', default=1000, type=int, help='Steps in integrator')
    parser.add_argument('--rng_seed', type=int, help='Seed for RNG', required=False,
                        default=random.randint(1000, 10000))
    parser.add_argument('--output_prefix',
                        default='sim', type=str,
                        help='Prefix for output files')
    parser.add_argument('--output_file',
                        default='trjout.h5', type=str,
                        help='Name of output trajectory file')
    parser.add_argument('--thermostat',
                        default='lv',
                        choices=('lv', 'vr'),
                        help='Thermostat to use, lv: Langevine, vr: Stochastic velocity rescale')
    parser.add_argument('--barostat', default='lv', choices=('lv', 'br'),
                        help='Barostat to use, lv: Langevine, br: Berendsen')
    parser.add_argument('--barostat_tau', default=5.0, type=float,
                        help='Tau parameter for Berendsen barostat')
    parser.add_argument('--barostat_mass', default=50.0, type=float,
                        help='Mass parameter for Langevin barostat')
    parser.add_argument('--barostat_gammaP', default=1.0, type=float,
                        help='gammaP parameter for Langevin barostat')
    parser.add_argument('--thermostat_gamma', type=float, default=0.5,
                        help='Thermostat coupling constant')
    parser.add_argument('--temperature', default=423.0, type=float, help='Temperature')
    parser.add_argument('--pressure', help='Pressure', type=float)
    parser.add_argument('--trj_collect', default=1000, type=int,
                        help='Collect trajectory every (step)')
    parser.add_argument('--energy_collect', default=1000, type=int,
                        help='Collect energy every (step)')
    parser.add_argument('--dt', default=0.001, type=float,
                        help='Integrator time step')
    parser.add_argument('--lj_cutoff', default=1.2, type=float,
                        help='Cutoff of atomistic non-bonded interactions')
    parser.add_argument('--cg_cutoff', default=1.4, type=float,
                        help='Cuoff of coarse-grained non-bonded interactions')
    parser.add_argument('--coulomb_epsilon1', default=1.0, type=float,
                        help='Epsilon_1 for coulomb interactions')
    parser.add_argument('--coulomb_epsilon2', default=78.0, type=float,
                        help='Epsilon_2 for coulomb interactions')
    parser.add_argument('--coulomb_kappa', default=0.0, type=float,
                        help='Kappa paramter for coulomb interactions')
    parser.add_argument('--coulomb_cutoff', default=0.9, type=float,
                        help='Cut-off for generalized reactive coulomb reactions')
    parser.add_argument('--table_groups', default='A,B',
                        help='Name of CG groups to read from tables')
    parser.add_argument('--initial_step', default=0,
                        help='Initial integrator step (useful for continue simulation',
                        type=int)
    parser.add_argument('--reactions', default=None,
                        help='Configuration file with chemical reactions')
    parser.add_argument('--debug', default=None, help='Turn on logging mechanism')
    parser.add_argument('--start_ar', default=0, type=int, help='When to start chemical reactions')
    parser.add_argument('--interactive', default=False, type=ast.literal_eval,
                        help='Run interactive mode')
    parser.add_argument('--store_species', default=False, type=ast.literal_eval,
                        help='Store particle types')
    parser.add_argument('--store_state', default=True, type=ast.literal_eval,
                        help='Store chemical state')
    parser.add_argument('--store_lambda', default=False, type=ast.literal_eval,
                        help='Store lambda parameter')

    return parser


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
    parser.add_argument('--dt_dyn', default=0.001, type=float, help='Integrator time step during backmapping')
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
    parser.add_argument('--cap_force', default=1000.0, type=float,
        help='Max force or 0.0 to switch it off cap-force')
    parser.add_argument('--two_phase', default=False, type=ast.literal_eval,
                        help='Two phase process, first bonded terms and then non-bonded')

    return parser
