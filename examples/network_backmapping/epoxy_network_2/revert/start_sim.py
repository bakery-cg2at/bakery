import argparse
import espressopp  # NOQA
import math  # NOQA
try:
    import MPI
except ImportError:
    from mpi4py import MPI
import time

import h5md_analysis as serial_h5md
import tools

kb = 0.0083144621  # GROMACS, kJ/molK

dt = 0.001
lj_cutoff = 1.2
cg_cutoff = 1.2
max_cutoff = max([lj_cutoff, cg_cutoff])

## Mostly you do not need to modify lines below.

def _args():
    parser = argparse.ArgumentParser('Backmapping simulation')
    parser.add_argument('--conf', required=True, help='Coordinate file')
    parser.add_argument('--top', required=True, help='Topology file')
    parser.add_argument('--node_grid', help='Node grid configuration')
    parser.add_argument('--skin', type=float, help='Skin value for VerletList')
    parser.add_argument('--resolution', default=0.0, type=float, help='Initial resolution')
    parser.add_argument('--res_rate', type=float, default=0.1, help='Alpha parameter')
    parser.add_argument('--eq', type=int, default=10000, help='CG run time')
    parser.add_argument('--thermostat_gamma', type=float, default=0.5)
    parser.add_argument('--int_step', default=1000, type=int, help='Number of steps in integration phase')
    parser.add_argument('--long', default=10000, type=int, help='Number of steps after backmapping')
    parser.add_argument('--rng_seed', default=12345, type=int, help='Seed for random number generator')
    parser.add_argument('--output_prefix', default='', type=str)
    parser.add_argument('--thermostat', default='lv', choices=('lv', 'vr'))
    parser.add_argument('--temperature', default=423.0, type=float)
    return parser.parse_args()


def main():  #NOQA
    args = _args()

    time0 = time.time()
    input_conf = espressopp.tools.gromacs.read(args.conf, args.top)

    N_atoms = len(input_conf.types)
    density = sum(input_conf.masses) / (input_conf.Lx * input_conf.Ly * input_conf.Lz)
    box = (input_conf.Lx, input_conf.Ly, input_conf.Lz)
    print('Setting up simulation...')
    print('Density = {}'.format(density))
    print('Box = {}'.format(box))

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    k_eq_step = int(args.eq/integrator_step)
    long_step = int(args.long/integrator_step)
    dynamic_res_time = int(int(1.0/args.res_rate)/integrator_step) if args.res_rate > 0.0 else 0
    sim_step = dynamic_res_time + k_eq_step + long_step
    end_dynamic_res_time = k_eq_step + dynamic_res_time
    if end_dynamic_res_time == k_eq_step:
        end_dynamic_res_time += 1

    if args.skin:
        skin = args.skin
    else:
        skin = 0.3

    print('Running with box {}'.format(box))
    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(args.rng_seed))

    system = tools.System(kb=kb)
    system.rng = espressopp.esutil.RNG(args.rng_seed)
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

    system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)
    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = dt

    part_prop, all_particles, adress_tuple = tools.genParticleList(
        input_conf, use_velocity=True, adress=True, use_charge=True)
    particle_ids = [x[0] for x in all_particles]
    print('Reads {} particles with properties {}'.format(len(all_particles), part_prop))

    system.storage.addParticles(all_particles, *part_prop)
    adress_fixed_list = espressopp.FixedTupleListAdress(system.storage)
    adress_fixed_list.addTuples(adress_tuple)
    system.storage.setFixedTuplesAdress(adress_fixed_list)

    print("Added tuples, decomposing now ...")
    system.storage.decompose()

# Exclude all bonded interaction from the lennard jones
    exclusionlist = input_conf.exclusions
    print('Excluded pairs from LJ interaction: {}'.format(len(exclusionlist)))

    verletlist = espressopp.VerletListAdress(
        system,
        cutoff=cg_cutoff,
        adrcut=lj_cutoff,
        dEx=0.0,
        dHy=max(box),
        adrCenter=[0.5*max(box), 0.5*max(box), 0.5*max(box)],
        exclusionlist=exclusionlist
        )

    vl_interaction, bondedinteractions, angleinteractions, dihedralinteractions, pairinteractions, cg_vl_interaction = (
        {}, {}, {}, {}, {}, {})
    vl_interaction = tools.setLennardJonesInteractions(
        system, input_conf, verletlist, lj_cutoff, input_conf.nonbond_params,
        ftpl=adress_fixed_list)
    cg_vl_interaction = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlist, cutoff=cg_cutoff, interaction=vl_interaction, ftpl=adress_fixed_list)
    #coulombinteraction = espressopp.tools.gromacs.setCoulombInteractions(
    #    system, verletlist, lj_cutoff, input_conf.atomtypeparams, epsilon1=1, epsilon2=80, kappa=0,
    #    adress=True, ftpl=adress_fixed_list)
    bondedinteractions = tools.setBondedInteractions(
        system, input_conf, ftpl=adress_fixed_list)
    #angleinteractions = tools.setAngleInteractions(
    #    system, input_conf, ftpl=adress_fixed_list)
    #dihedralinteractions = tools.setDihedralInteractions(
    #    system, input_conf, ftpl=adress_fixed_list)
    #pairinteractions = tools.setPairInteractions(
    #    system, input_conf, lj_cutoff, ftpl=adress_fixed_list)

    print('='*10)
    print('Bonds: {}'.format(sum(len(x) for x in input_conf.bondtypes.values())))
    print('Angles: {}'.format(sum(len(x) for x in input_conf.angletypes.values())))
    print('Dihedrals: {}'.format(sum(len(x) for x in input_conf.dihedraltypes.values())))
    print('Pairs: {}'.format(sum(len(x) for x in input_conf.pairtypes.values())))
    print('='*10)

    dynamic_res = espressopp.integrator.DynamicResolution(
        system,
        verletlist,
        adress_fixed_list,
        args.res_rate)
    integrator.addExtension(dynamic_res)
    dynamic_res.active = False

# Define the thermostat
    if args.temperature:
        temperature = args.temperature*kb
    print('Temperature: {}, gamma: {}'.format(temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    if args.thermostat == 'lv':
        thermostat = espressopp.integrator.LangevinThermostat(system)
        thermostat.temperature = temperature
        thermostat.adress = True
        thermostat.gamma = args.thermostat_gamma
    elif args.thermostat == 'vr':
        thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
        thermostat.temperature = temperature
        thermostat.coupling = args.thermostat_gamma
    integrator.addExtension(thermostat)

    print("Added tuples, decomposing now ...")
    system.storage.decompose()
    print('AdResS decomposition')
    espressopp.tools.AdressDecomp(system, integrator)
# Write the warmup configuration

# --- below is the simulation --- #
# warmup
    print('Trajectory saved to: {}but.h5'.format(args.output_prefix))
    #h5dump = serial_h5md.DumpH5MD(
    #    '{}melf.h5'.format(args.output_prefix),
    #    system,
    #    integrator,
    #    edges=list(box),
    #    particle_ids=particle_ids,
    #    unfolded=True,
    #    save_vel=False
    #    )

    #h5dump.dump()
    #h5dump.analyse()
    #h5dump.flush()

    print('Energy saved to: {}energy_{}_.csv'.format(args.output_prefix, args.res_rate))
    system_analysis = espressopp.analysis.SystemAnalysis(
        system,
        integrator,
        '{}energy_{}_.csv'.format(args.output_prefix, args.res_rate))
    system_analysis.add_observable('res', espressopp.analysis.Resolution(system, dynamic_res))
    system_analysis.add_observable('lj', espressopp.analysis.PotentialEnergy(system, vl_interaction))
    for (lb, cross), interHarmonic in bondedinteractions.iteritems():
        system_analysis.add_observable('bond_%d%s' % (lb, '_cross' if cross else ''),
                                       espressopp.analysis.PotentialEnergy(system, interHarmonic), False)
    for (lb, cross), interAngularHarmonic in angleinteractions.iteritems():
        system_analysis.add_observable('angle_%d%s' % (lb, '_cross' if cross else ''),
                                       espressopp.analysis.PotentialEnergy(system, interAngularHarmonic), False)
    for (lb, cross), interDihHarmonic in dihedralinteractions.iteritems():
       system_analysis.add_observable('dihedral_%d%s' % (lb, '_cross' if cross else ''),
                                      espressopp.analysis.PotentialEnergy(system, interDihHarmonic), False)
    for (lb, cross), interaction14 in pairinteractions.iteritems():
        system_analysis.add_observable('lj-14_%d%s' % (lb, '_cross' if cross else ''),
                                       espressopp.analysis.PotentialEnergy(system, interaction14), False)
    #system_analysis.add_observable('coul', espressopp.analysis.PotentialEnergy(system, coulombinteraction), True)

    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, 1)
    integrator.addExtension(ext_analysis)

    system_analysis.dump()

    traj_gro = espressopp.io.DumpGROAdress(
        system, adress_fixed_list, integrator, filename='traj.gro', unfolded=True, append=True)

    print('Simulation for steps: %d' % (sim_step*integrator_step))
    print('Dynamic resolution, rate={}'.format(args.res_rate))
    print('CG equilibration for {}'.format(k_eq_step*integrator_step))
    print('Measuring energy with higher resolution for {}'.format(
        (end_dynamic_res_time-k_eq_step)*integrator_step))
    print('Long run for {}'.format(long_step*integrator_step))

    system_analysis.info()
    for k in range(sim_step):
        if k == k_eq_step:
            print('End of CG simulation. Start dynamic resolution')
            dynamic_res.active = True
            ext_analysis.interval = 1
        if k == end_dynamic_res_time:
            print('End of dynamic resolution, change energy measuring accuracy to 500')
            ext_analysis.interval = 1
        integrator.run(integrator_step)
        # total_velocity.reset()
        #if k % 100 == 0:
        #    h5dump.dump()
        #    h5dump.flush()
        traj_gro.dump()
        system_analysis.info()
    #h5dump.dump()
    #h5dump.close()
    system_analysis.info()
    traj_gro.dump()

    print('finished!')
    print('total time: {}'.format(time.time() - time0))
    espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())


if __name__ == '__main__':
    main()
