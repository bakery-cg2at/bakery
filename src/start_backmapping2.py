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
import random
import os
import time

import files_io
import tools_sim as tools
import gromacs_topology

import tools_backmapping

from app_args import _args_backmapping as _args

# GROMACS units, kJ/mol K
kb = 0.0083144621

# Storage options
simulation_author = os.environ.get('USER', 'xxx')
simulation_email = 'xxx@xxx.xxx'


# Mostly you do not need to modify lines below.

def main():  # NOQA
    h5md_group = 'atoms'
    print h5md_group
    time0 = time.time()
    args = _args().parse_args()
    _args().save_to_file('{}params.out'.format(args.output_prefix), args)

    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([args.lj_cutoff, args.cg_cutoff])

    print('Welcome in bakery!\n')

    print('Reading hybrid topology and coordinate file')
    input_conf = gromacs_topology.read(args.conf, args.top)
    input_gro_conf = files_io.GROFile(args.conf)
    input_gro_conf.read()

    box = (input_conf.Lx, input_conf.Ly, input_conf.Lz)
    print('\nSetting up simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    k_eq_step = int(args.eq / integrator_step)
    long_step = int(args.long / integrator_step)
    dynamic_res_time = int(int(1.0 / args.alpha) / integrator_step) + 1 if args.alpha > 0.0 else 0

    if args.skin:
        skin = args.skin
    else:
        skin = 0.16

    rng_seed = args.rng_seed
    if args.rng_seed == -1:
        rng_seed = random.randint(1, 10000)

    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(rng_seed))
    print('Time step: {}'.format(args.dt))
    print('LJ cutoff: {}'.format(lj_cutoff))
    print('CG cutoff: {}'.format(cg_cutoff))
    print('Boltzmann constant = {}'.format(kb))

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG(rng_seed)

    part_prop, all_particles, adress_tuple = tools.genParticleList(
        input_conf, use_velocity=True, adress=True, use_charge=True)
    print('Reads {} particles with properties {}'.format(len(all_particles), part_prop))

    # Make output from AT particles.
    at_gro_conf = files_io.GROFile.copy(input_gro_conf, [x for p in adress_tuple for x in p[1:]], renumber=True)
    gro_whole = files_io.GROFile.copy(input_gro_conf, [x for p in adress_tuple for x in p], renumber=True)

    # Generate initial velocities, only for CG particles.
    particle_list = []
    index_adrat = part_prop.index('adrat')
    if 'v' not in part_prop:
        print('Generating velocities from Maxwell-Boltzmann distribution for T={}'.format(
            args.temperature))
        part_prop.append('v')
        cg_particles = [x for x in all_particles if x.adrat == 0]
        vx, vy, vz = espressopp.tools.velocities.gaussian(
            args.temperature, len(cg_particles), [x.mass for x in cg_particles],
            kb=kb)
        cg_id = 0
        last_vel = (0.0, 0.0, 0.0)
        for p in all_particles:
            t = list(p)
            if p.adrat == 0:
                last_vel = (vx[cg_id], vy[cg_id], vz[cg_id])
                cg_id += 1
            del t[index_adrat]
            t.append(espressopp.Real3D(last_vel))
            particle_list.append(t)
    else:
        for p in all_particles:
            t = list(p)
            del t[index_adrat]
            particle_list.append(t)

    del part_prop[index_adrat]

    print('Running with box {}'.format(box))
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    if args.node_grid:
        nodeGrid = map(int, args.node_grid.split(','))
    else:
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    print('Number of nodes {}, node-grid: {}'.format(
        MPI.COMM_WORLD.size, nodeGrid))
    if args.cell_grid:
        cellGrid = map(int, args.cell_grid.split(','))
    else:
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, max_cutoff, skin)
    print('Cell grid: {}'.format(cellGrid))

    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    system.storage.addParticles(map(tuple, particle_list), *part_prop)

    vs_list = espressopp.FixedVSList(system.storage)
    vs_list.addTuples(adress_tuple)

    at_particle_group = espressopp.ParticleGroup(system.storage)
    at_particle_ids = set()
    cg_particle_ids = set()
    for a in adress_tuple:
        cg_particle_ids.add(a[0])
        for at in a[1:]:
            at_particle_group.add(at)
            at_particle_ids.add(at)

    integrator = espressopp.integrator.VelocityVerletHybrid(system, vs_list)
    integrator.dt = args.dt

    system.storage.decompose()

    print('Prepared:')
    print('Bonds: {}'.format(sum(len(x) for x in input_conf.bondtypes.values())))
    print('Angles: {}'.format(sum(len(x) for x in input_conf.angletypes.values())))
    print('Dihedrals: {}'.format(sum(len(x) for x in input_conf.dihedraltypes.values())))
    print('Pairs: {}'.format(sum(len(x) for x in input_conf.pairtypes.values())))

    print('Setting dynamic resolution')

    dynamic_res = espressopp.integrator.DynamicResolution(
        system,
        vs_list,
        args.alpha)
    integrator.addExtension(dynamic_res)
    dynamic_res.active = False
    dynamic_res.resolution = args.initial_resolution

    # Define interactions.
    verletlistAT = None
    verletlistCG = None
    if args.two_phase:
        verletlistCG = tools_backmapping.setupFirstPhase(
            system, args, input_conf, at_particle_ids, cg_particle_ids)
    else:
        verletlistAT, verletlistCG = tools_backmapping.setupSinglePhase(
            system, args, input_conf, at_particle_ids, cg_particle_ids)

    # Define the thermostat
    if args.temperature:
        temperature = args.temperature * kb
    print('Temperature: {}, gamma: {}'.format(args.temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    thermostat = espressopp.integrator.LangevinThermostatOnGroup(system, at_particle_group)
    thermostat.temperature = temperature
    thermostat.gamma = args.thermostat_gamma
    integrator.addExtension(thermostat)

    print("Added tuples, decomposing now ...")

    output_file = 'trjout.h5'
    h5file = '{}_{}_{}_{}'.format(
        args.output_prefix,
        rng_seed, args.alpha,
        output_file)
    print('Trajectory saved to: {}'.format(h5file))
    traj_file = espressopp.io.DumpH5MD(
        system, h5file,
        group_name=h5md_group,
        static_box=False,
        author=simulation_author,
        email=simulation_email,
        store_lambda=True,
        store_species=True)
    traj_file.set_parameters({
        'temperature': temperature,
        'thermostat': args.thermostat,
        'skin': skin,
        'rng_seed': rng_seed,
        'lj_cutoff': lj_cutoff,
        'cg_cutoff': cg_cutoff,
        'vl_cutoff': max_cutoff,
        'integrator_step': integrator_step,
        'dt': args.dt
    })

    # Sets energy storage.
    energy_file = '{}energy_{}_{}.csv'.format(args.output_prefix, args.alpha, args.rng_seed)
    print('Energy saved to: {}'.format(energy_file))
    system_analysis = espressopp.analysis.SystemMonitor(
        system,
        integrator,
        espressopp.analysis.SystemMonitorOutputCSV(energy_file))
    temp_comp = espressopp.analysis.Temperature(system)
    system_analysis.add_observable('T', temp_comp)
    system_analysis.add_observable('Ekin', espressopp.analysis.KineticEnergy(system, temp_comp))
    system_analysis.add_observable(
        'res', espressopp.analysis.Resolution(system, dynamic_res))

    system_info_filter = args.system_info_filter.split(',') if args.system_info_filter else None
    for label, interaction in sorted(system.getAllInteractions().items()):
        show_in_system_info = True
        if system_info_filter is not None:
            show_in_system_info = False
            for v in system_info_filter:
                if v in label:
                    show_in_system_info = True
                    break
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction),
            show_in_system_info)

    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, args.energy_collect)
    integrator.addExtension(ext_analysis)

    system_analysis.dump()

    k_trj_collect = int(math.ceil(float(args.trj_collect) / integrator_step))

    print('Dynamic resolution, rate={}'.format(args.alpha))
    print('CG equilibration for {}'.format(k_eq_step * integrator_step))
    print('Collect trajectory every {} step'.format(k_trj_collect * integrator_step))
    print('Collect energy every {} step'.format(args.energy_collect))
    print('Atomistic long run for {}'.format(long_step * integrator_step))

    traj_file.dump(integrator.step, integrator.step * args.dt)

    has_capforce = False
    if args.cap_force > 0.0:
        has_capforce = True
        print('Define maximum cap-force during the backmapping')
        cap_force = espressopp.integrator.CapForce(system, args.cap_force)

    system.storage.decompose()

    system_analysis.info()

    ############# SIMULATION #####################
    global_int_step = 0
    # First run the eq phase.
    for k in range(k_eq_step):
        integrator.run(integrator_step)
        system_analysis.info()
        traj_file.dump(global_int_step * integrator_step, global_int_step * integrator_step * args.dt)
        global_int_step += 1
    else:
        global_int_step += 1
        system_analysis.info()
        traj_file.dump(global_int_step*integrator_step, global_int_step*integrator_step*args.dt)

    traj_file.flush()
    traj_file.close()
    # Now run backmapping.
    dynamic_res.active = True
    print('Change time-step to {}'.format(args.dt_dyn))
    integrator.dt = args.dt_dyn
    ext_analysis.interval = args.energy_collect_bck
    if has_capforce:
        integrator.addExtension(cap_force)
    print('End of CG simulation. Start dynamic resolution, dt={}'.format(
        args.dt_dyn))
    if args.two_phase:
        # Run first phase, only bonded terms and non-bonded CG term is enabled.
        for k in range(dynamic_res_time):
            integrator.run(integrator_step)
            system_analysis.info()
            global_int_step += 1

        print('First phase finished, switch time-step to {}'.format(args.dt))
        integrator.dt = args.dt
        print('Simulate for {} steps'.format(10*integrator_step))
        for k in range(10):
            integrator.run(integrator_step)
            system_analysis.info()
            global_int_step += 1

        confout_aa = '{}confout_aa_{}_{}_phase_one.gro'.format(args.output_prefix, args.alpha, args.rng_seed)
        at_gro_conf.update_positions(system)
        at_gro_conf.write(confout_aa, force=True)
        print('Atomistic configuration write to: {}'.format(confout_aa))

	########## SECOND PHASE ################
        # Change interactions.
        print('Second phase, switch on non-bonded interactions, time-step: {}'.format(args.dt_dyn))
        verletlistCG.disconnect()
        verletlistAT, verletlistCG = tools_backmapping.setupSecondPhase(
            system, args, input_conf, at_particle_ids, cg_particle_ids)
        # Reset dynamic res, start again.
        if args.alpha2 is not None:
            print('Change dynamic resolution alpha: '.format(args.alpha2))
            dynamic_res.rate = args.alpha2
        dynamic_res.active = True
        dynamic_res.resolution = args.initial_resolution

        # Reset system analysis.
        ext_analysis.disconnect()

        energy_file = '{}energy_{}_{}_two.csv'.format(args.output_prefix, args.alpha, args.rng_seed)
        print('Energy saved to: {}'.format(energy_file))
        system_analysis2 = espressopp.analysis.SystemMonitor(
            system,
            integrator,
            espressopp.analysis.SystemMonitorOutputCSV(energy_file))
        system_analysis2.copy_state(system_analysis)
        system_analysis2.add_observable('T', temp_comp)
        system_analysis2.add_observable('Ekin', espressopp.analysis.KineticEnergy(system, temp_comp))
        system_analysis2.add_observable(
            'res', espressopp.analysis.Resolution(system, dynamic_res))

        for label, interaction in sorted(system.getAllInteractions().items()):
            show_in_system_info = True
            if system_info_filter is not None:
                show_in_system_info = False
                for v in system_info_filter:
                    if v in label:
                        show_in_system_info = True
                        break
            system_analysis2.add_observable(
                label, espressopp.analysis.PotentialEnergy(system, interaction), show_in_system_info)

        ext_analysis2 = espressopp.integrator.ExtAnalyze(
            system_analysis2, args.energy_collect_bck)
        integrator.addExtension(ext_analysis2)

        # Simulation
        for k in range(dynamic_res_time+10):
            integrator.run(integrator_step)
            system_analysis2.info()
            global_int_step += 1
        else:
            system_analysis2.info()

        if has_capforce:
            print('Switch off cap-force')
            cap_force.disconnect()
    else:
        for k in range(dynamic_res_time):
            integrator.run(integrator_step)
            system_analysis.info()
            traj_file.dump(global_int_step * integrator_step,
                           global_int_step * integrator_step * args.dt)
            global_int_step += 1
        else:
            global_int_step += 1
            system_analysis.info()
            traj_file.dump(global_int_step * integrator_step,
                           global_int_step * integrator_step * args.dt)
            traj_file.close()

    confout_aa = '{}confout_aa_{}_{}_phase_two.gro'.format(args.output_prefix, args.alpha, args.rng_seed)
    at_gro_conf.update_positions(system)
    at_gro_conf.write(confout_aa, force=True)

    print('Atomistic configuration write to: {}'.format(confout_aa))

    ############ Now run normal AT simulation.############
    print('End of dynamic resolution, change energy measuring accuracy to {}'.format(
        args.energy_collect))
    print('Set back time-step to: {}'.format(args.dt))
    ext_analysis.interval = args.energy_collect
    if args.two_phase:
        ext_analysis2.interval = args.energy_collect
    integrator.dt = args.dt
    print('Running for {} steps'.format(long_step * integrator_step))
    for k in range(long_step):
        integrator.run(integrator_step)
    for k in range(long_step):
        global_int_step += 1
        if args.two_phase:
            system_analysis2.info()
        else:
            system_analysis.info()
    else:
        global_int_step += 1
        if args.two_phase:
            system_analysis2.info()
        else:
            system_analysis.info()
        traj_file.dump(global_int_step*integrator_step, global_int_step*integrator_step*args.dt)
        traj_file.dump(global_int_step * integrator_step, global_int_step * integrator_step * args.dt)
        traj_file.close()

    confout_aa = '{}confout_aa_{}_{}.gro'.format(args.output_prefix, args.alpha, args.rng_seed)
    at_gro_conf.update_positions(system)
    at_gro_conf.write(confout_aa, force=True)
    print('Atomistic configuration write to: {}'.format(confout_aa))

    print('Finished!')
    print('Total time: {}'.format(time.time() - time0))
    espressopp.tools.analyse.final_info(system, integrator, verletlistAT, time0, time.time())


if __name__ == '__main__':
    main()
