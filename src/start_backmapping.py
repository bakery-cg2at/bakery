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
import logging

import espressopp  # NOQA
import math  # NOQA
from mpi4py import MPI
import random
import os
import time

from . import files_io
from . import tools_sim as tools
from . import gromacs_topology

from . import tools_backmapping
from . import tools as general_tools

from .app_args import _args_backmapping as _args

from .logger import logger

# GROMACS units, kJ/mol K
kb = 0.0083144621

# Storage options
simulation_author = os.environ.get('USER', 'xxx')
simulation_email = 'xxx@xxx.xxx'

# Mostly you do not need to modify lines below.

def main(args):  # NOQA
    h5md_group = 'atoms'
    time0 = time.time()

    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([args.lj_cutoff, args.cg_cutoff])

    print('Welcome in bakery!\n')

    print('Reading hybrid topology and coordinate file')

    generate_exclusions = False #args.exclusion_list is None or not os.path.exists(args.exclusion_list)      

    input_conf = gromacs_topology.read(args.top, doRegularExcl=generate_exclusions)
    input_gro_conf = files_io.GROFile(args.conf)
    input_gro_conf.read()

    if not generate_exclusions:
        exclusion_file = open(args.exclusion_list, 'r')
        exclusions = [list(map(int, x.split())) for x in exclusion_file.readlines()]
        print(('Read exclusion list from {} (total: {})'.format(args.exclusion_list, len(exclusions))))
        input_conf = input_conf._replace(exclusions=exclusions)
    else:
        exclusion_list_file = 'exclusion_{}.list'.format(args.top.split('.')[0])
        with open(exclusion_list_file, 'w') as fel:
            for p in input_conf.exclusions:
                fel.write('{} {}\n'.format(*p))
        print(('Save exclusion list: {} ({})'.format(exclusion_list_file, len(input_conf.exclusions))))

    box = input_gro_conf.box
    print('\nSetting up simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    if args.trj_collect > 0:
        integrator_step = min([integrator_step, args.trj_collect])
    k_eq_step = int(args.eq / integrator_step)
    long_step = int(args.long / integrator_step)
    dynamic_res_time = 0
    if args.alpha > 0.0:
        dynamic_res_time = int(int(1.0 / args.alpha) / integrator_step) + 2
        if args.nonuniform_lambda:
            dynamic_res_time += int(10000/integrator_step)
            print(('Running nonuniform lambda, extended running time by {} steps'.format(10000)))

    if args.skin:
        skin = args.skin
    else:
        skin = 0.16

    rng_seed = args.rng_seed
    if args.rng_seed == -1:
        rng_seed = random.randint(1, 10000)
        args.rng_seed = rng_seed

    random.seed(rng_seed)

    _args().save_to_file('{}_{}_params.out'.format(args.output_prefix, rng_seed), args)

    print(('Skin: {}'.format(skin)))
    print(('RNG Seed: {}'.format(rng_seed)))
    print(('Time step: {}'.format(args.dt)))
    print(('LJ cutoff: {}'.format(lj_cutoff)))
    print(('CG cutoff: {}'.format(cg_cutoff)))
    print(('Boltzmann constant = {}'.format(kb)))

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG(rng_seed)

    part_prop, all_particles, adress_tuple = tools.genParticleList(
        input_conf, input_gro_conf, adress=True, use_charge=True)
    print(('Reads {} particles with properties {}'.format(len(all_particles), part_prop)))

    if input_conf.charges:
        print(('Total charge: {}'.format(sum(input_conf.charges))))

    # Make output from AT particles.
    at_gro_conf = files_io.GROFile.copy(input_gro_conf, [x for p in adress_tuple for x in p[1:]], renumber=True)
    gro_whole = files_io.GROFile.copy(input_gro_conf, [x for p in adress_tuple for x in p], renumber=True)

    # Generate initial velocities, only for CG particles, AT particles will get the CG particle velocity.
    particle_list = []
    index_adrat = part_prop.index('adrat')
    print(('Generating velocities from Maxwell-Boltzmann distribution for T={}'.format(
        args.temperature)))
    part_prop.append('v')
    cg_particles = [x for x in all_particles if x.adrat == 0]
    vx, vy, vz = espressopp.tools.velocities.gaussian(
        args.temperature,
        len(cg_particles),
        [x.mass for x in cg_particles],
        kb=kb)
    cg_id = 0
    last_vel = (0.0, 0.0, 0.0)
    last_lambda = 0.0
    last_res_id = -1
    index_lambda = part_prop.index('lambda_adr')
    for p in all_particles:
        t = list(p)
        if p.adrat == 0:  # this is CG particle
            last_vel = (vx[cg_id], vy[cg_id], vz[cg_id])
            cg_id += 1
        if args.nonuniform_lambda and p.res_id != last_res_id:
            last_res_id = p.res_id
            last_lambda = -1.0 * random.uniform(0.0, 10000 * args.alpha)
        if args.nonuniform_lambda:
            t[index_lambda] = last_lambda
        del t[index_adrat]
        t.append(espressopp.Real3D(last_vel))
        particle_list.append(t)

    del part_prop[index_adrat]

    print(('Running with box {}'.format(box)))
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    if args.node_grid:
        nodeGrid = list(map(int, args.node_grid.split(',')))
    else:
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    print(('Number of nodes {}, node-grid: {}'.format(
        MPI.COMM_WORLD.size, nodeGrid)))
    if args.cell_grid:
        cellGrid = list(map(int, args.cell_grid.split(',')))
    else:
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, max_cutoff, skin)
    print(('Cell grid: {}'.format(cellGrid)))

    hook_after_init = lambda *_, **__: True
    hook_setup_interactions = lambda *_, **__: True
    if args.hooks and os.path.exists(args.hooks):
        print(('Found {}'.format(args.hooks)))
        l = {}
        exec(compile(open(args.hooks, "rb").read(), args.hooks, 'exec'), globals(), l)
        hook_after_init = l.get('hook_after_init', hook_after_init)
        hook_setup_interactions = l.get('hook_setup_interactions', hook_setup_interactions)

    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    hook_after_init(system, particle_list, adress_tuple, input_conf, part_prop)

    system.storage.addParticles(list(map(tuple, particle_list)), *part_prop)
    system.storage.decompose()

    vs_list = espressopp.FixedVSList(system.storage)
    vs_list.addTuples(adress_tuple)

    at_particle_group = espressopp.ParticleGroup(system.storage)
    cg_particle_group = espressopp.ParticleGroup(system.storage)
    particle_groups = {'at': at_particle_group, 'cg': cg_particle_group}
    at_particle_ids = set()
    cg_particle_ids = set()
    for a in adress_tuple:
        cg_particle_ids.add(a[0])
        cg_particle_group.add(a[0])
        for at in a[1:]:
            at_particle_group.add(at)
            at_particle_ids.add(at)

    integrator = espressopp.integrator.VelocityVerletHybrid(system, vs_list)
    integrator.dt = args.dt

    system.integrator = integrator

    system.storage.decompose()

    print('Prepared:')
    print(('Bonds: {}'.format(sum(len(x) for x in list(input_conf.bondtypes.values())))))
    print(('Angles: {}'.format(sum(len(x) for x in list(input_conf.angletypes.values())))))
    print(('Dihedrals: {}'.format(sum(len(x) for x in list(input_conf.dihedraltypes.values())))))
    print(('Pairs: {}'.format(sum(len(x) for x in list(input_conf.pairtypes.values())))))
    print(('CG particles: {}'.format(len(cg_particle_ids))))
    print(('AT particles: {}'.format(len(at_particle_ids))))

    print('Setting dynamic resolution')
    dynamic_res = espressopp.integrator.DynamicResolution(
        system,
        vs_list,
        args.alpha)
    integrator.addExtension(dynamic_res)
    dynamic_res.active = False
    dynamic_res.resolution = args.initial_resolution

    if args.table_groups is None:
        table_groups = []
    else:
        table_groups = args.table_groups.split(',')
        print(('Using table groups: {}'.format(table_groups)))

    # Define interactions.
    verletlistAT, verletlistCG = tools_backmapping.setupSinglePhase(
        system, args, input_conf, at_particle_ids, cg_particle_ids, table_groups=table_groups)

    hook_setup_interactions(system, input_conf, verletlistAT, verletlistCG)

    print(('Number of interactions: {}'.format(system.getNumberOfInteractions())))

    # Define the thermostat
    temperature = args.temperature * kb
    print(('Temperature: {} ({}), gamma: {}'.format(args.temperature, temperature, args.thermostat_gamma)))
    print(('Thermostat: {}'.format(args.thermostat)))
    if args.thermostat == 'lv':
        if args.thermostat_whole:
            print('Enable thermostat on all particles, not only atomistic')
            thermostat = espressopp.integrator.LangevinThermostat(system)
        elif args.thermostat_cg:
            print('Enable thermostat only on CG particles')
            thermostat = espressopp.integrator.LangevinThermostatOnGroup(system, cg_particle_group)
        else:
            thermostat = espressopp.integrator.LangevinThermostatOnGroup(system, at_particle_group)
        thermostat.temperature = temperature
        thermostat.gamma = args.thermostat_gamma
    elif args.thermostat == 'vr':
        thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
        thermostat.temperature = temperature
        thermostat.coupling = args.thermostat_gamma
    integrator.addExtension(thermostat)

    print("Added tuples, decomposing now ...")

    output_file = 'trjout.h5'
    h5file = '{}_{}_{}_{}'.format(
        args.output_prefix,
        rng_seed, args.alpha,
        output_file)
    print(('Trajectory saved to: {}'.format(h5file)))
    traj_file = espressopp.io.DumpH5MD(
        system, h5file,
        group_name=h5md_group,
        static_box=True,
        author=simulation_author,
        email=simulation_email,
        store_lambda=True,
        store_species=True,
        store_force=args.store_force,
        store_state=args.store_state,
        is_single_prec=True,
        chunk_size=256)
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

    ext_analysis, system_analysis = tools.setSystemAnalysis(
        system, integrator, args, args.energy_collect, '_first', dynamic_res, particle_groups)
    system_analysis.dump()

    k_trj_collect = int(math.ceil(float(args.trj_collect) / integrator_step))
    k_trj_flush = 25 if 25 < 10*k_trj_collect else 10*k_trj_collect
    if k_trj_collect == 0:
        k_trj_flush = 0

    print(('Dynamic resolution, rate={}'.format(args.alpha)))
    print(('CG equilibration for {}'.format(k_eq_step * integrator_step)))
    print(('Collect trajectory every {} step'.format(k_trj_collect * integrator_step)))
    print(('Flush trajectory every {} step'.format(k_trj_flush * integrator_step)))
    print(('Collect energy every {} step'.format(args.energy_collect)))
    print(('Atomistic long run for {}'.format(long_step * integrator_step)))

    system.storage.decompose()

    if args.gro_collect > 0:
        gro_collect_filename = '{}confout_dump_{}_{}.gro'.format(
            args.output_prefix, args.alpha, rng_seed)
        dump_conf_gro = espressopp.io.DumpGRO(system, integrator, filename=gro_collect_filename, append=True)
        ext_dump_conf_gro = espressopp.integrator.ExtAnalyze(
            dump_conf_gro, args.gro_collect)
        integrator.addExtension(ext_dump_conf_gro)
        print(('Store .gro file {}'.format(gro_collect_filename)))

    if args.remove_com > 0:
        print(('Removes total velocity of the system every {} steps'.format(args.remove_com)))
        total_velocity = espressopp.analysis.CMVelocity(system)
        ext_remove_com = espressopp.integrator.ExtAnalyze(total_velocity, args.remove_com)
        integrator.addExtension(ext_remove_com)

    gro_whole.update_positions(system)
    gro_whole.write(
        '{}confout_full_{}_{}_before.gro'.format(args.output_prefix, args.alpha, args.rng_seed), force=True)

    print(('Number of particles: {}'.format(len(particle_list))))

    time_sim0 = time.time()

    ############# SIMULATION: EQUILIBRATION PHASE #####################
    system_analysis.dump()
    global_int_step = 0
    for k in range(k_eq_step):
        system_analysis.info()
        if k_trj_collect > 0 and k % k_trj_collect == 0:
            traj_file.dump(global_int_step * integrator_step, global_int_step * integrator_step * args.dt)
        if k_trj_flush > 0 and k % k_trj_flush == 0:
            traj_file.flush()  # Write HDF5 to disk.
        integrator.run(integrator_step)
        global_int_step += 1

    time_cg = time.time() - time_sim0
    system_analysis.dump()

    ######### Now run backmapping.  #######################
    time_sim0 = time.time()
    has_capforce = False
    if args.cap_force and not args.cap_force_lj:
        has_capforce = True
        print(('Define maximum cap-force during the backmapping (max: {})'.format(args.cap_force)))
        cap_force = espressopp.integrator.CapForce(system, args.cap_force)

    print('Activating dynamic resolution changer')
    dynamic_res.active = True

    print(('Change time-step to {}'.format(args.dt_dyn)))
    integrator.dt = args.dt_dyn
    if has_capforce:
        thermostat.disconnect()
        integrator.addExtension(cap_force)
        thermostat.connect()
    print(('End of CG simulation. Start dynamic resolution, dt={}'.format(
        args.dt_dyn)))
    two_phase = args.two_phase or args.second_phase_em

    if two_phase:
        ext_analysis.disconnect()
        verletlistCG.disconnect()
        verletlistAT.disconnect()
        verletlistCG = tools_backmapping.setupFirstPhase(
            system, args, input_conf, at_particle_ids, cg_particle_ids)

        ext_analysis2, system_analysis2 = tools.setSystemAnalysis(
            system,
            integrator,
            args,
            args.energy_collect_bck,
            '_one',
            dynamic_res,
            particle_groups)

        # Run first phase, only bonded terms and non-bonded CG term is enabled.
        for k in range(dynamic_res_time):
            if k_trj_collect > 0 and k % k_trj_collect == 0:
                traj_file.dump(global_int_step * integrator_step, global_int_step * integrator_step * args.dt)
            if k_trj_flush > 0 and k % k_trj_flush == 0:
                traj_file.flush()  # Write HDF5 to disk.
            system_analysis2.info()
            integrator.run(integrator_step)
            global_int_step += 1

        confout_aa = '{}confout_aa_{}_{}_phase_one.gro'.format(args.output_prefix, args.alpha, rng_seed)
        at_gro_conf.update_positions(system)
        at_gro_conf.write(confout_aa, force=True)
        gro_whole.update_positions(system)
        gro_whole.write(
            '{}confout_full_{}_{}_phase_one.gro'.format(args.output_prefix, args.alpha, args.rng_seed), force=True)
        print(('Atomistic configuration write to: {}'.format(confout_aa)))

        ########## SECOND PHASE ################
        # Change interactions.
        print(('Second phase, switch on non-bonded interactions, time-step: {}'.format(args.dt_dyn)))
        verletlistCG.disconnect()
        verletlistAT, verletlistCG = tools_backmapping.setupSecondPhase(
            system, args, input_conf, at_particle_ids, cg_particle_ids)
        # Reset dynamic res, start again.
        if args.alpha2 is not None:
            print(('Change dynamic resolution alpha: {}'.format(args.alpha2)))
            dynamic_res.rate = args.alpha2
            dynamic_res_time = int(int(1.0 / args.alpha2) / integrator_step) + 1 if args.alpha2 > 0.0 else 0
        dynamic_res.active = True
        dynamic_res.resolution = args.initial_resolution

        # Reset system analysis.
        ext_analysis2.disconnect()

        ext_analysis3, system_analysis3 = tools.setSystemAnalysis(
            system, integrator, args, args.energy_collect_bck, '_two', dynamic_res, particle_groups)

        if args.second_phase_em:
            minimize_energy = espressopp.integrator.MinimizeEnergy(system, 0.0001, 10.0, 0.001 * input_conf.Lx)
            while not minimize_energy.run(100, True):
                pass
        else:
            # Simulation
            for k in range(dynamic_res_time):
                if k_trj_collect > 0 and k % k_trj_collect == 0:
                    traj_file.dump(global_int_step * integrator_step, global_int_step * integrator_step * args.dt)
                if k_trj_flush > 0 and k % k_trj_flush == 0:
                    traj_file.flush()  # Write HDF5 to disk.
                system_analysis3.info()
                integrator.run(integrator_step)
                global_int_step += 1
    else:
        # Single phase backmapping
        ext_analysis.interval = args.energy_collect_bck
        print('Running a single-phase backmapping.')
        for k in range(dynamic_res_time):
            if k_trj_collect > 0 and k % k_trj_collect == 0:
                traj_file.dump(global_int_step * integrator_step, global_int_step * integrator_step * args.dt)
            if k_trj_flush > 0 and k % k_trj_flush == 0 and k > 0:
                traj_file.flush()  # Write HDF5 to disk.
            system_analysis.info()
            integrator.run(integrator_step)
            global_int_step += 1

    # After backmapping, switch off dynamic resolution
    print('Disconnect dynamic_res')
    dynamic_res.active = False

    time_bck = time.time() - time_sim0

    gro_whole.update_positions(system)
    gro_whole.write(
        '{}confout_full_{}_{}_phase_two.gro'.format(args.output_prefix, args.alpha, rng_seed), force=True)
    confout_aa = '{}confout_aa_{}_{}_phase_two.gro'.format(args.output_prefix, args.alpha, rng_seed)
    at_gro_conf.update_positions(system)
    at_gro_conf.write(confout_aa, force=True)

    print(('Atomistic configuration write to: {}'.format(confout_aa)))

    ############ Now run normal AT simulation.############
    print(('End of dynamic resolution, change energy measuring accuracy to {}'.format(
        args.energy_collect)))
    print(('Set back time-step to: {}'.format(args.dt)))
    time_sim0 = time.time()
    ext_analysis.interval = args.energy_collect
    if two_phase:
        ext_analysis3.interval = args.energy_collect
    else:
        ext_analysis.interval = args.energy_collect

    if has_capforce:
        if args.cap_force_ramp is None:
            cap_force.disconnect()
            print('Cap-force switched off')
        else:
            cap_force.ramp = args.cap_force_ramp
            print(('Cap-force switched gradually, decrease of {}'.format(cap_force.ramp)))

    integrator.dt = args.dt
    print(('Running for {} steps'.format(long_step * integrator_step)))
    for k in range(long_step):
        if k_trj_collect > 0 and k % k_trj_collect == 0:
            traj_file.dump(global_int_step * integrator_step, global_int_step * integrator_step * args.dt)
        if k_trj_flush > 0 and k % k_trj_flush == 0:
            traj_file.flush()  # Write HDF5 to disk.
        if two_phase:
            system_analysis3.info()
        else:
            system_analysis.info()
        integrator.run(integrator_step)
        global_int_step += 1

    if args.em > 0:
        if has_capforce:
            cap_force.disconnect()
        print('Runninng basic energy minimization')
        if two_phase:
            system_analysis3.info()
        else:
            system_analysis.info()
        minimize_energy = espressopp.integrator.MinimizeEnergy(system, args.em_gamma, args.em_ftol, args.em_max_d * input_gro_conf.box[0], True)
        minimize_energy.run(args.em, True)
        print('Energy information:')
        if two_phase:
            system_analysis3.info()
        else:
            system_analysis.info()

    time_at = time.time() - time_sim0

    ## Save benchmark data
    if os.path.exists('{}benchmark.dat'.format(args.output_prefix)):
        benchmark_file = open('{}benchmark.dat'.format(args.output_prefix), 'a')
    else:
        benchmark_file = open('{}benchmark.dat'.format(args.output_prefix), 'w')
        benchmark_file.write('N_at\tN_cg\tCPUs\talpha\ttime_cg\ttime_bck\ttime_at\n')

    benchmark_file.write('{Nat}\t{Ncg}\t{CPUs}\t{alpha}\t{time_cg}\t{time_bck}\t{time_at}\n'.format(
        Nat=len(at_particle_ids), Ncg=len(cg_particle_ids), CPUs=MPI.COMM_WORLD.size, alpha=args.alpha,
        time_cg=time_cg,time_bck=time_bck,time_at=time_at))
    benchmark_file.close()
    ## End save benchmark data

    gro_whole.update_positions(system)
    gro_whole.write(
        '{}confout_final_full_{}_{}.gro'.format(args.output_prefix, args.alpha, rng_seed), force=True)

    confout_aa = '{}confout_final_aa_{}_{}.gro'.format(args.output_prefix, args.alpha, rng_seed)
    at_gro_conf.update_positions(system)
    at_gro_conf.write(confout_aa, force=True)
    print(('Final atomistic configuration write to: {}'.format(confout_aa)))
    print(('Final hybrid configuration write to: {}'.format(
        '{}confout_final_full_{}_{}.gro'.format(args.output_prefix, args.alpha, rng_seed))))

    # Write atomistic topology
    hyb_top = files_io.GROMACSTopologyFile(args.top)
    hyb_top.read()
    at_topology = general_tools.get_atomistic_topology(
        hyb_top,
        virtual_atomtypes=[
            v['atnum'] for v in list(input_conf.atomtypeparams.values()) if v['particletype'] == 'V'])
    topol_aa = '{}topol_final_aa_{}_{}.top'.format(args.output_prefix, args.alpha, rng_seed)
    at_topology.write(topol_aa)
    print(('Final AA topology: {}'.format(topol_aa)))

    traj_file.close()

    print('Finished!')
    print(('Total time: {}'.format(time.time() - time0)))
    espressopp.tools.analyse.final_info(system, integrator, verletlistAT, time0, time.time())


if __name__ == '__main__':
    args = _args().parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.ERROR)

    if args.debug:
        import ipdb
        with ipdb.launch_ipdb_on_exception():
            main(args)
    else:
        main(args)
