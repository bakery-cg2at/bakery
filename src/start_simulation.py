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
import shutil

import numpy as np
import tools_sim as tools
import tools as general_tools
import gromacs_topology

from app_args import _args_md

# GROMACS units, kJ/mol K
kb = 0.0083144621

h5md_group = 'atoms'

__doc__ = 'Run GROMACS-like simulation'

# Do not to modify lines below.


def sort_trajectory(trj, ids):
    print('Sorting trajectory')
    idd = [
        x[1] for x in sorted([(p_id, col_id) for col_id, p_id in enumerate(ids)],
                             key=lambda y: (True, y[0]) if y[0] == -1 else (False, y[0]))
    ]
    return trj[idd]


def main():  #NOQA
    args = _args_md().parse_args()

    _args_md().save_to_file('{}params.out'.format(args.output_prefix), args)

    if args.debug:
        for s in args.debug.split(','):
            print('Activating logger {}'.format(s))
            logging.getLogger(s.strip()).setLevel(logging.DEBUG)

    table_groups = map(str.strip, args.table_groups.split(','))
    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([lj_cutoff, cg_cutoff])
    dt = args.dt

    time0 = time.time()
    input_conf = gromacs_topology.read(args.conf, args.top)

    box = (input_conf.Lx, input_conf.Ly, input_conf.Lz)
    print('Setup simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    sim_step = args.run / integrator_step

    if args.skin:
        skin = args.skin

    rng_seed = args.rng_seed
    if not args.rng_seed:
        rng_seed = random.randint(10, 1000000)

    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(rng_seed))

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

    if args.coord:
        import h5py
        print("Reading coordinates from {}".format(args.coord))
        h5coord = h5py.File(args.coord)
        pos = h5coord['/particles/{}/position/value'.format(h5md_group)][args.coord_frame]
        try:
            ids = h5coord['/particles/{}/id/value'.format(h5md_group)][args.coord_frame]
            has_ids = True
        except:
            has_ids = False
        pos = h5coord['/particles/{}/position/value'.format(h5md_group)][args.coord_frame]
        if has_ids:
            pos = sort_trajectory(pos, ids)

        try:
            species = h5coord['/particles/{}/species/value'.format(h5md_group)][args.coord_frame]
        except:
            species = h5coord['/particles/{}/species'.format(h5md_group)]
        if has_ids:
            species = sort_trajectory(species, ids)
        try:
            box = np.array(
                h5coord['/particles/{}/box/edges/value'.format(h5md_group)][args.coord_frame])
        except:
            box = np.array(h5coord['/particles/{}/box/edges'.format(h5md_group)])
        valid_species = {x[1] for x in all_particles}
        ppid = 0
        for pid, p in enumerate(pos):
            if species[pid] in valid_species:
                particle_list[ppid][2] = espressopp.Real3D(p)
                ppid += 1
        h5coord.close()

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
    vl_interaction = gromacs_topology.setLennardJonesInteractions(
        system, input_conf.defaults, input_conf.atomtypeparams,
        verletlist, lj_cutoff, input_conf.nonbond_params, table_groups=table_groups)
    gromacs_topology.setTabulatedInteractions(
        system, input_conf.atomtypeparams, vl=verletlist,
        cutoff=cg_cutoff, interaction=vl_interaction, table_groups=table_groups)
    bondedinteractions = gromacs_topology.setBondedInteractions(
        system, input_conf.bondtypes, input_conf.bondtypeparams)
    angleinteractions = gromacs_topology.setAngleInteractions(
        system, input_conf.angletypes, input_conf.angletypeparams)
    dihedralinteractions = gromacs_topology.setDihedralInteractions(
        system, input_conf.dihedraltypes, input_conf.dihedraltypeparams)
    pairinteractions = gromacs_topology.setPairInteractions(
        system, input_conf.pairtypes, input_conf.pairtypeparams, lj_cutoff)
    coulomb_interaction = gromacs_topology.setCoulombInteractions(
        system, verletlist, args.coulomb_cutoff, input_conf.atomtypeparams,
        epsilon1=args.coulomb_epsilon1,
        epsilon2=args.coulomb_epsilon2,
        kappa=args.coulomb_kappa)

    if coulomb_interaction: 
        system.addInteraction(coulomb_interaction, 'coulomb')

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
    energy_file = '{}_energy_{}.csv'.format(args.output_prefix, rng_seed)
    print('Energy saved to: {}'.format(energy_file))
    system_analysis = espressopp.analysis.SystemMonitor(
        system,
        integrator,
        espressopp.analysis.SystemMonitorOutputCSV(energy_file))
    temp_comp = espressopp.analysis.Temperature(system)
    system_analysis.add_observable('T', temp_comp)
    system_analysis.add_observable(
        'Ekin', espressopp.analysis.KineticEnergy(
            system, temp_comp))
    for label, interaction in sorted(system.getAllInteractions().items()):
        print('System analysis: adding {}'.format(label))
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction))
    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, args.energy_collect)
    integrator.addExtension(ext_analysis)
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
    dump_gro = espressopp.io.DumpGRO(
        system, integrator, filename=output_gro_file,
        unfolded=True, append=False)
    dump_gro.dump()
    print('Wrote end configuration to: {}'.format(output_gro_file))

    print('finished!')
    print('total time: {}'.format(time.time()-time0))
    espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())


if __name__ == '__main__':
    main()
