"""
Copyright (C) 2016,2019 Jakub Krajniak <jkrajniak@gmail.com>

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


def _args_analyze():
    parser = general_tools.MyArgParser(description='Analyze GROMACS topology',
                                       fromfile_prefix_chars='@')
    parser.add_argument('--top', '--topology', required=True, help='Topology file')
    
    return parser

def _args_backmapping():
    parser = general_tools.MyArgParser(
        description='Starts MD simulation and run backmapping',
        fromfile_prefix_chars='@')
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--verbose', default=False, action='store_true')
    parser.add_argument('--conf', required=True, help='Coordinate file')
    parser.add_argument('--top', '--topology', required=True, help='Topology file',
                        dest='top')
    parser.add_argument('--node_grid', help='Node grid configuration')
    parser.add_argument('--cell_grid', help='Cell grid configuration')
    parser.add_argument('--skin', type=float, help='Skin value for VerletList')
    parser.add_argument('--initial_resolution',
                        default=0.0, type=float, help='Initial resolution')
    parser.add_argument('--alpha',
                        type=float, default=0.1, help='Alpha parameter')
    parser.add_argument('--alpha2',
                        type=float, default=None, help='Alpha parameter (second phase)')
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
    parser.add_argument('--thermostat_whole',
                        default=False, type=ast.literal_eval,
                        help='Thermalize all particles, not only AT')
    parser.add_argument('--thermostat_cg',
                        default=False, type=ast.literal_eval,
                        help='Thermalize all particles, only CG')
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
    parser.add_argument('--coulomb_epsilon2', default=78.0, type=float,
                        help='Epsilon 2 parameter for Coulomb reaction field potential')
    parser.add_argument('--coulomb_kappa', default=0.0, type=float,
                        help='Kappa paramter for coulomb interactions')
    parser.add_argument('--coulomb_cutoff', default=0.9, type=float,
                        help='Cutoff for coulomb interaction')
    parser.add_argument('--system_info_filter', default=None, type=str)
    parser.add_argument('--energy_collect', default=1000, type=int,
                        help='Collect energy every (step)')
    parser.add_argument('--energy_collect_bck', default=1, type=int,
                        help='Collect energy every (step) during backmapping')
    parser.add_argument('--trj_collect', default=1000, type=int,
                        help='Collect trajectory every (step)')
    parser.add_argument('--gro_collect', default=0, type=int, help='If set then collect trajcectory in .gro')
    parser.add_argument('--cap_force', default=None, type=float,
        help='Max force or 0.0 to switch it off cap-force')
    parser.add_argument('--cap_force_ramp', default=None, type=float,
                        help='Gradually switch off cap-force')
    parser.add_argument('--cap_force_lj', default=500000.0, type=float, help='Max force only for non-bonded terms')
    parser.add_argument('--two_phase', default=False, type=ast.literal_eval,
                        help='Two phase process, first bonded terms and then non-bonded')
    parser.add_argument('--second_phase_em', default=False, type=ast.literal_eval,
                        help='Second phase with minimize energy')
    parser.add_argument('--exclusion_list', default='exclusion_hyb_topol.list', required=True,
                        help='The exclusion list')
    parser.add_argument('--remove_com', type=int, default=0,
                        help='Resets the total velocity of the system every n-th steps')
    parser.add_argument('--store_species', default=False, type=ast.literal_eval,
                        help='Store particle types')
    parser.add_argument('--store_force', default=False, type=ast.literal_eval,
                        help='Store particle force')
    parser.add_argument('--store_state', default=False, type=ast.literal_eval,
                        help='Store chemical state')
    parser.add_argument('--store_lambda', default=True, type=ast.literal_eval,
                        help='Store lambda parameter')
    parser.add_argument('--table_groups', required=True,
                        help='Name of CG groups to read from tables')
    parser.add_argument('--nonuniform_lambda', default=False, type=ast.literal_eval,
                        help='Distribute initial lambda non-uniformly in the box')
    parser.add_argument('--save_interactions', default=False, type=ast.literal_eval)
    parser.add_argument('--hooks', default=None)

    parser.add_argument('--disable_h5md', default=False, type=ast.literal_eval)

    disable_interactions = parser.add_argument_group('Disabled interactions')
    disable_interactions.add_argument('--disable_angles', default=False, type=ast.literal_eval)
    disable_interactions.add_argument('--disable_dihedrals', default=False, type=ast.literal_eval)
    
    em_group = parser.add_argument_group('Energy minimization')
    em_group.add_argument('--em', help='Maximum number of steps to perform in EM', type=int, default=0)
    em_group.add_argument('--em_gamma', help='Gamma parameter for force damping', type=float, default=0.0001)
    em_group.add_argument('--em_ftol', help='Force tolerance', type=float, default=10.0)
    em_group.add_argument('--em_max_d', help='Max displacement x box dimension', default=0.001, type=float)



    return parser
