#!/usr/bin/env python
"""
Copyright (C) 2015-2016 Jakub Krajniak <jkrajniak@gmail.com>

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

import espressopp  # NOQA
import math  # NOQA



try:
    import MPI
except ImportError:
    from mpi4py import MPI
import time
import logging
import random

import files_io
import gromacs_topology

from app_args import _args_md

import tools_sim as tools

# GROMACS units, kJ/mol K
kb = 0.0083144621

h5md_group = 'atoms'

__doc__ = 'Run GROMACS-like simulation'

# Basically, you do not need to modify lines below.

def main():  #NOQA
    args = _args_md().parse_args()

    if args.debug:
        for s in args.debug.split(','):
            print('Activating logger {}'.format(s))
            logging.getLogger(s.strip()).setLevel(logging.DEBUG)

    if args.table_groups:
        table_groups = map(str.strip, args.table_groups.split(','))
    else:
        table_groups = []
    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([lj_cutoff, cg_cutoff])
    dt = args.dt

    time0 = time.time()

    generate_exclusions = args.exclusion_list is None

    input_conf = gromacs_topology.read(args.conf, args.top, doRegularExcl=generate_exclusions)
    input_conf_gro = files_io.GROFile(args.conf)
    input_conf_gro.read()

    if generate_exclusions:
        # Save exclusions for future runs.
        exclusion_file = open('exclusions_{}.list'.format(args.top.split('.')[0]), 'w')
        for e in input_conf.exclusions:
            exclusion_file.write('{} {}\n'.format(*e))
        exclusion_file.close()
        print('Exclusion list written down')
    else:
        exclusion_file = open(args.exclusion_list, 'r')
        exclusions = [map(int, x.split()) for x in exclusion_file.readlines()]
        print('Read exclusion list from {} (total: {})'.format(args.exclusion_list, len(exclusions)))
        input_conf = input_conf._replace(exclusions=exclusions)

    box = (input_conf.Lx, input_conf.Ly, input_conf.Lz)
    print('Setup simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    sim_step = args.run / integrator_step

    if args.skin:
        skin = args.skin

    rng_seed = args.rng_seed
    if args.rng_seed == -1:
        rng_seed = random.randint(1, 10000)
        args.rng_seed = rng_seed

    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(rng_seed))

    _args_md().save_to_file('{}_{}_params.out'.format(args.output_prefix, rng_seed), args)

    part_prop, all_particles = gromacs_topology.genParticleList(
        input_conf, use_velocity=True, use_charge=True)
    print('Reads {} particles with properties {}'.format(len(all_particles), part_prop))

    particle_list = []
    if 'v' not in part_prop:
        print('Generating velocities from Maxwell-Boltzmann distribution for T={}'.format(
            args.temperature))
        part_prop.append('v')
        vx, vy, vz = espressopp.tools.velocities.gaussian(
            args.temperature, len(all_particles), [x.mass for x in all_particles],
            kb=kb)
        ppid = 0
        for p in all_particles:
            t = list(p)
            t.append(espressopp.Real3D(vx[ppid], vy[ppid], vz[ppid]))
            particle_list.append(t)
    else:
        particle_list = map(list, all_particles)

    density = sum(input_conf.masses)*1.6605402 / (box[0] * box[1] * box[2])
    print('Density: {} kg/m^3'.format(density))
    print('Box: {} nm'.format(box))

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG(rng_seed)
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    if args.node_grid:
        nodeGrid = map(int, args.node_grid.split(','))
    else:
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    print('Number of nodes {}, node-grid: {}'.format(
        MPI.COMM_WORLD.size, nodeGrid))
    cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, max_cutoff, skin)

    print('Cell grid: {}'.format(cellGrid))

    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = dt
    system.integrator = integrator

    system.storage.addParticles(particle_list, *part_prop)
    system.storage.decompose()

# In the case of butane is very easy to do
    print('Excluded pairs from LJ interaction: {}'.format(len(input_conf.exclusions)))

# Exclude all bonded interaction from the lennard jones
    verletlist = espressopp.VerletList(
        system,
        cutoff=max_cutoff,
        exclusionlist=input_conf.exclusions
        )

# define the potential, interaction_id = 0
    lj_interaction = tools.setLennardJonesInteractions(
        system, input_conf, verletlist, lj_cutoff,
        input_conf.nonbond_params,
        interaction=espressopp.interaction.VerletListLennardJones(verletlist),
        table_groups=table_groups)
    if lj_interaction is not None:
        system.addInteraction(lj_interaction, 'lj')

    if args.coulomb_cutoff > 0.0:
        coulomb_interaction = gromacs_topology.setCoulombInteractions(
            system, verletlist, args.coulomb_cutoff, input_conf.atomtypeparams,
            epsilon1=args.coulomb_epsilon1,
            epsilon2=args.coulomb_epsilon2, kappa=args.coulomb_kappa,
            interaction=espressopp.interaction.VerletListReactionFieldGeneralized(
                verletlist))
        if coulomb_interaction is not None:
            system.addInteraction(coulomb_interaction, 'coulomb')

    tools.setBondedInteractions(
        system, input_conf)
    tools.setAngleInteractions(
        system, input_conf)
    tools.setDihedralInteractions(
        system, input_conf)
    tools.setPairInteractions(
        system, input_conf, args.lj_cutoff)
    tab_cg = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlist,
        cutoff=cg_cutoff,
        interaction=espressopp.interaction.VerletListTabulated(
            verletlist,
        ), table_groups=table_groups)
    if tab_cg is not None:
        system.addInteraction(tab_cg, 'tab-cg')

    print('Bonds: {}'.format(sum(len(x) for x in input_conf.bondtypes.values())))
    print('Angles: {}'.format(sum(len(x) for x in input_conf.angletypes.values())))
    print('Dihedrals: {}'.format(sum(len(x) for x in input_conf.dihedraltypes.values())))
    print('Pairs: {}'.format(sum(len(x) for x in input_conf.pairtypes.values())))

# Define the thermostat
    temperature = args.temperature*kb
    print('Temperature: {} ({}), gamma: {}'.format(temperature, temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    if args.thermostat == 'lv':
        thermostat = espressopp.integrator.LangevinThermostat(system)
        thermostat.temperature = temperature
        thermostat.gamma = args.thermostat_gamma
    elif args.thermostat == 'vr':
        thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
        thermostat.temperature = temperature
        thermostat.coupling = args.thermostat_gamma
    else:
        raise Exception('Wrong thermostat keyword: `{}`'.format(args.thermostat))
    integrator.addExtension(thermostat)

    pressure_comp = espressopp.analysis.Pressure(system)
    if args.pressure:
        pressure = args.pressure * 0.060221374  # convert from bars to gromacs units kj/mol/nm^3
        if args.barostat == 'lv':
            print('Barostat: Langevin with P={}, gamma={}, mass={}'.format(
                pressure, 0.5, pow(10, 4)))
            barostat = espressopp.integrator.LangevinBarostat(system, system.rng, temperature)
            barostat.gammaP = args.barostat_gammaP
            barostat.mass = args.barostat_mass
            barostat.pressure = pressure
        elif args.barostat == 'br':
            print('Barostat: Berendsen with P={} and tau={}'.format(pressure, 0.5))
            barostat = espressopp.integrator.BerendsenBarostat(system, pressure_comp)
            barostat.tau = args.barostat_tau
            barostat.pressure = pressure
        integrator.addExtension(barostat)

    print("Decomposing now ...")
    system.storage.decompose()

    # Observe tuple lists
    args.alpha = 0.0
    ext_analysis, system_analysis = tools.setSystemAnalysis(
        system, integrator, args, args.energy_collect, '')
    print('Configured system analysis')

    h5md_output_file = '{}_{}_{}'.format(args.output_prefix, rng_seed, args.output_file)
    print('Save trajectory to: {}'.format(h5md_output_file))
    traj_file = espressopp.io.DumpH5MD(
        system, h5md_output_file,
        group_name=h5md_group,
        static_box=False,
        author='Jakub Krajniak',
        email='jkrajniak@gmail.com',
        store_species=args.store_species,
        store_state=args.store_state,
        store_lambda=args.store_lambda)

    gro_whole = files_io.GROFile.copy(input_conf_gro)
    file_name = '{}_{}_confout_time_0.gro'.format(args.output_prefix, rng_seed)
    print('Write coordinates before to: {}'.format(file_name))
    gro_whole.update_positions(system)
    gro_whole.write(file_name, force=True)

    traj_file.set_parameters({'temperature': args.temperature})

    print('Reset total velocity')
    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    integrator.step = args.initial_step

    k_trj_collect = int(math.ceil(float(args.trj_collect) / integrator_step))
    k_energy_collect = int(math.ceil(float(args.energy_collect) / integrator_step))

    print('Running simulation for {} steps'.format(sim_step*integrator_step))
    print('Collect trajectory every {} step'.format(k_trj_collect*integrator_step))
    print('Collect energy every {} step'.format(k_energy_collect*integrator_step))

    if args.interactive:
        import IPython
        IPython.embed()

    system_analysis.dump()
    # Save interactions

    if args.em > 0:
        print('Runninng basic energy minimization')
        system_analysis.info()
        minimize_energy = espressopp.integrator.MinimizeEnergy(system, args.em_gamma, args.em_ftol, args.em_max_d * input_conf.Lx)
        minimize_energy.run(args.em, True)
        print('Energy information:')
        system_analysis.dump()
        system_analysis.info()
    else:
        for k in range(sim_step):
            if k_energy_collect > 0 and k % k_energy_collect == 0:
                system_analysis.info()
            if k_trj_collect > 0 and k % k_trj_collect == 0:
                int_step = args.initial_step + k*integrator_step
                traj_file.dump(int_step, int_step*dt)
            if k_trj_collect > 0 and k % 100 == 0:
                traj_file.flush()
            integrator.run(integrator_step)
        else:
            traj_file.dump(sim_step*integrator_step, sim_step*integrator_step*dt)
            traj_file.close()

    # Saves output file.
    output_gro_file = '{}_{}_confout.gro'.format(args.output_prefix, rng_seed)
    gro_whole.update_positions(system)
    gro_whole.write(output_gro_file, force=True)
    print('Wrote end configuration to: {}'.format(output_gro_file))
    print('Save interaction database...')
    #tools.saveInteractions(system, '{}_{}_interactions.pck'.format(args.output_prefix, rng_seed))

    print('finished!')
    print('total time: {}'.format(time.time()-time0))
    #espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())


if __name__ == '__main__':
    main()
