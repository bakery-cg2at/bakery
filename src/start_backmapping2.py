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

import tools_sim as tools
import gromacs_topology

from app_args import _args_backmapping as _args

# GROMACS units, kJ/mol K
kb = 0.0083144621

# Storage options
simulation_author = os.environ.get('USER', 'xxx')
simulation_email = 'xxx@xxx.xxx'

# Mostly you do not need to modify lines below.



def main():  #NOQA
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

    box = (input_conf.Lx, input_conf.Ly, input_conf.Lz)
    print('\nSetting up simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    k_eq_step = int(args.eq/integrator_step)
    long_step = int(args.long/integrator_step)
    dynamic_res_time = int(int(1.0/args.alpha)/integrator_step) if args.alpha > 0.0 else 0
    sim_step = dynamic_res_time + k_eq_step + long_step
    end_dynamic_res_time = k_eq_step + dynamic_res_time
    if end_dynamic_res_time == k_eq_step:
        end_dynamic_res_time += 1

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

    if args.coord:
        import h5py
        print("Reading coordinates from {}".format(args.coord))
        h5coord = h5py.File(args.coord)
        h5md_group = args.coord_h5md_group
        pos = h5coord['/particles/{}/position/value'.format(h5md_group)][-1]
        try:
            species = h5coord['/particles/{}/species/value'.format(h5md_group)][-1]
        except:
            species = h5coord['/particles/{}/species'.format(h5md_group)]
        try:
            box = list(h5coord['/particles/{}/box/edges/value'.format(h5md_group)][-1])
        except:
            box = list(h5coord['/particles/{}/box/edges'.format(h5md_group)])
        valid_species = {x[1] for x in all_particles}
        ppid = 0
        for pid, p in enumerate(pos):
            if species[pid] in valid_species:
                particle_list[ppid][2] = espressopp.Real3D(p)
                ppid += 1
        h5coord.close()

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
    for a in adress_tuple:
        for at in a[1:]:
            at_particle_group.add(at)

    integrator = espressopp.integrator.VelocityVerletHybrid(system, vs_list)
    integrator.dt = args.dt

    system.storage.decompose()

    try:
        os.remove('details.gro')
        print('removed old details.gro')
    except OSError:
        pass
    dump_gro_details = espressopp.io.DumpGRO(
        system,
        integrator,
        filename='details.gro',
        unfolded=True,
        append=True)
    ext_dump_details = espressopp.integrator.ExtAnalyze(dump_gro_details, 100)
    dump_gro_details.dump()
    integrator.addExtension(ext_dump_details)

# Exclude all bonded interaction from the lennard jones
    exclusionlist = input_conf.exclusions
    print('Excluded pairs from LJ interaction: {}'.format(len(exclusionlist)))
    verletlistAT = espressopp.VerletListHybridAT(
        system, cutoff=lj_cutoff, exclusionlist=exclusionlist)

    verletlistCG = espressopp.VerletListHybridCG(
        system, cutoff=cg_cutoff, exclusionlist=exclusionlist)

    lj_interaction = espressopp.interaction.VerletListHybridLennardJones(
        verletlistAT, False)
    tools.setLennardJonesInteractions(
        system, input_conf, verletlistAT, lj_cutoff,
        input_conf.nonbond_params,
        interaction=lj_interaction)
    coulomb_interaction = espressopp.interaction.VerletListHybridReactionFieldGeneralized(
        verletlistAT, False)
    gromacs_topology.setCoulombInteractions(
        system, verletlistAT, lj_cutoff, input_conf.atomtypeparams,
        epsilon1=args.coulomb_epsilon1,
        epsilon2=args.coulomb_epsilon2, kappa=args.coulomb_kappa,
        interaction=coulomb_interaction)
    tools.setBondedInteractions(
        system, input_conf)
    tools.setAngleInteractions(
        system, input_conf)
    tools.setDihedralInteractions(
        system, input_conf)
    tools.setPairInteractions(
        system, input_conf, lj_cutoff)
    tab_cg = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlistCG,
        cutoff=cg_cutoff,
        interaction=espressopp.interaction.VerletListHybridTabulated(
            verletlistCG, True
        ))
    system.addInteraction(lj_interaction, 'xyz-lj')
    system.addInteraction(coulomb_interaction, 'xyz-coulomb')
    system.addInteraction(tab_cg, 'xyz-cg')
    print lj_interaction
    lj_interaction.scale_factor = 0.0
    print coulomb_interaction
    coulomb_interaction.scale_factor = 0.0

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

# Define the thermostat
    if args.temperature:
        temperature = args.temperature*kb
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

    for label, interaction in sorted(system.getAllInteractions().items()):
        print('System analysis: adding {}'.format(label))
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction))

    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, args.energy_collect)
    integrator.addExtension(ext_analysis)

    system_analysis.dump()

    k_trj_collect = int(math.ceil(float(args.trj_collect) / integrator_step))

    print('Dynamic resolution, rate={}'.format(args.alpha))
    print('CG equilibration for {}'.format(k_eq_step*integrator_step))
    print('Measuring energy with higher resolution for {}'.format(
        (end_dynamic_res_time-k_eq_step)*integrator_step))
    print('Collect trajectory every {} step'.format(k_trj_collect*integrator_step))
    print('Atomistic long run for {}'.format(long_step*integrator_step))
    print('Running for {} steps...'.format(sim_step*integrator_step))

    traj_file.dump(integrator.step, integrator.step*args.dt)

    cap_force = espressopp.integrator.CapForce(system, 1000.)
    integrator.addExtension(cap_force)

    system.storage.decompose()

    system_analysis.info()
    for k in range(sim_step):
        if k == k_eq_step:
            print('End of CG simulation. Start dynamic resolution.')
            dynamic_res.active = True
            ext_analysis.interval = args.energy_collect_bck
        if k == end_dynamic_res_time:
            print('End of dynamic resolution, change energy measuring accuracy to {}'.format(
                args.energy_collect))
            ext_analysis.interval = args.energy_collect
            integrator.dt = 0.001
            print('Set LJ and Coulombic interaction with full scale')
            lj_interaction.scale_factor = 1.0
            coulomb_interaction.scale_factor = 1.0

        integrator.run(integrator_step)

        # total_velocity.reset()
        if k_trj_collect > 0 and k % k_trj_collect == 0:
            traj_file.dump(k*integrator_step, k*integrator_step*args.dt)
        if k_trj_collect > 0 and k % 100 == 0:
            traj_file.flush()
        system_analysis.info()
    else:
        system_analysis.dump()
        system_analysis.info()
        traj_file.dump(sim_step*integrator_step, sim_step*integrator_step*args.dt)
        traj_file.close()

    confout_aa = '{}confout_aa_{}_{}.gro'.format(args.output_prefix, args.alpha, args.rng_seed)
    dump_gro = espressopp.io.DumpGRO(
        system,
        integrator,
        filename=confout_aa,
        unfolded=True,
        append=False)
    dump_gro.dump()

    print('Finished!')
    print('Atomistic configuration write to: {}'.format(confout_aa))
    print('Total time: {}'.format(time.time() - time0))
    espressopp.tools.analyse.final_info(system, integrator, verletlistAT, time0, time.time())


if __name__ == '__main__':
    main()
