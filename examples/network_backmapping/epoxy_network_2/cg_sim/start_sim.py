import argparse
import espressopp  # NOQA
import math  # NOQA
try:
    import MPI
except ImportError:
    from mpi4py import MPI
import time

import tools

kb = 0.0083144621  # GROMACS, kJ/molK

dt = 0.001
max_cutoff = 1.4
lj_cutoff = 1.4
cg_cutoff = 1.4
table_groups = ['A', 'B', 'C', 'D', 'E', 'F']

h5md_group = 'atoms'
filename = 'but.h5'

# Do not to modyfied lines below.

def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--conf', required=True, help='Input .gro coordinate file')
    parser.add_argument('--top', required=True, help='Topology file')
    parser.add_argument('--node_grid')
    parser.add_argument('--skin', type=float, default=0.16)
    parser.add_argument('--coord', help='Input coordinate h5md file')
    parser.add_argument('--coord_frame', default=-1, type=int)
    parser.add_argument('--run', type=int, default=10000)
    parser.add_argument('--int_step', default=1000, type=int, help='Steps in integrator')
    parser.add_argument('--rng_seed', default=12345, type=int, help='Seed for RNG')
    parser.add_argument('--output_prefix', default='', type=str)
    parser.add_argument('--output_file', default='but.h5', type=str)
    parser.add_argument('--thermostat', default='lv', choices=('lv', 'vr'))
    parser.add_argument('--barostat', default='lv', choices=('lv', 'br'))
    parser.add_argument('--thermostat_gamma', type=float, default=0.5)
    parser.add_argument('--temperature', default=423.0, type=float, help='Temperature (K)')
    parser.add_argument('--pressure', help='Pressure (bar)', type=float)
    parser.add_argument('--trj_collect', default=1000, type=int,
                        help='Collect trajectory every (step)')
    parser.add_argument('--energy_collect', default=1000, type=int,
                        help='Collect energy every (step)')
    parser.add_argument('--settings', help='Settings file')
    return parser.parse_args()


def process_config(args):
    config = {}
    if args.settings:
        execfile(args.settings, {}, config)
    config.update(vars(args))
    return config


def main():  #NOQA
    args = _args()

    time0 = time.time()
    input_conf = espressopp.tools.gromacs.read(args.conf, args.top)

    N_atoms = len(input_conf.types)
    density = sum(input_conf.masses)*1.6605402 / (input_conf.Lx * input_conf.Ly * input_conf.Lz)
    box = (input_conf.Lx, input_conf.Ly, input_conf.Lz)
    print('Setup simulation...')
    print('Density: {} kg/m^3'.format(density))
    print('Box: {} nm'.format(box))

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    sim_step = args.run / integrator_step

    if args.skin:
        skin = args.skin

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

    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = dt
    system.integrator = integrator

    part_prop, all_particles = espressopp.tools.gromacs.genParticleList(input_conf, use_velocity=True)
    print('Reads {} particles with properties {}'.format(len(all_particles), part_prop))

    if args.coord:
        import h5py
        print("Reading coordinates from {}".format(args.coord))
        h5coord = h5py.File(args.coord)
        pos = h5coord['/particles/{}/position/value'.format(h5md_group)][-1]
        species = h5coord['/particles/{}/species'.format(h5md_group)]
        valid_species = {x[1] for x in all_particles}
        ppid = 0
        for pid, p in enumerate(pos):
            if species[pid] in valid_species:
                all_particles[ppid][2] = espressopp.Real3D(p)
                ppid += 1
        h5coord.close()

    particle_ids = [x[0] for x in all_particles]

    system.storage.addParticles(all_particles, *part_prop)
    system.storage.decompose()

# In the case of butane is very easy to do
    exclusionlist = input_conf.exclusions

    print('Excluded pairs from LJ interaction: {}'.format(len(exclusionlist)))

# Exclude all bonded interaction from the lennard jones
    verletlist = espressopp.VerletList(
        system,
        cutoff=max_cutoff,
        exclusionlist=exclusionlist
        )

# define the potential, interaction_id = 0
    vl_interaction = espressopp.tools.gromacs.setLennardJonesInteractions(
        system, input_conf.defaults,input_conf.atomtypeparams,
        verletlist, lj_cutoff, input_conf.nonbond_params, table_groups=table_groups)
    cg_vl_interaction = espressopp.tools.gromacs.setTabulatedInteractions(
        system, input_conf.atomtypeparams, vl=verletlist,
        cutoff=cg_cutoff, interaction=vl_interaction, table_groups=table_groups)
    bondedinteractions = espressopp.tools.gromacs.setBondedInteractions(
        system, input_conf.bondtypes, input_conf.bondtypeparams)
    angleinteractions = espressopp.tools.gromacs.setAngleInteractions(
        system, input_conf.angletypes, input_conf.angletypeparams)
    dihedralinteractions = espressopp.tools.gromacs.setDihedralInteractions(
        system, input_conf.dihedraltypes, input_conf.dihedraltypeparams)
    pairinteractions = espressopp.tools.gromacs.setPairInteractions(
        system, input_conf.pairtypes, input_conf.pairtypeparams, lj_cutoff)

    print('Bonds: {}'.format(sum(len(x) for x in input_conf.bondtypes.values())))
    print('Angles: {}'.format(sum(len(x) for x in input_conf.angletypes.values())))
    print('Dihedrals: {}'.format(sum(len(x) for x in input_conf.dihedraltypes.values())))
    print('Pairs: {}'.format(sum(len(x) for x in input_conf.pairtypes.values())))

# Define the thermostat
    temperature = args.temperature*kb
    print('Temperature: {}, gamma: {}'.format(temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    if args.thermostat == 'lv':
        thermostat = espressopp.integrator.LangevinThermostat(system)
        thermostat.temperature = temperature
        thermostat.gamma = args.thermostat_gamma
    elif args.thermostat == 'vr':
        thermostat = espressopp.integrator.StochasticVelocityRescaling(system)
        thermostat.temperature = temperature
        thermostat.coupling = args.thermostat_gamma
    integrator.addExtension(thermostat)

    pressure_comp = espressopp.analysis.Pressure(system)
    if args.pressure:
        pressure = args.pressure * 0.060221374  # convert from bars to gromacs units kj/mol/nm^3
        if args.barostat == 'lv':
            print('Barostat: Langevin with P={}, gamma={}, mass={}'.format(pressure, 0.5, pow(10, 4)))
            barostat = espressopp.integrator.LangevinBarostat(system, system.rng, temperature)
            barostat.gammaP = 0.5
            barostat.mass = pow(10, 4)
            barostat.pressure = pressure
        elif args.barostat == 'br':
            print('Barostat: Berendsen with P={} and tau={}'.format(pressure, 0.5))
            barostat = espressopp.integrator.BerendsenBarostat(system, pressure_comp)
            barostat.tau = 50.0
            barostat.pressure = pressure
        integrator.addExtension(barostat)

    print("Decomposing now ...")
    system.storage.decompose()

    print('Energy saved to: {}energy.csv'.format(args.output_prefix))
    system_analysis = espressopp.analysis.SystemAnalysis(
        system,
        integrator,
        '{}energy.csv'.format(args.output_prefix))
    system_analysis.add_observable('lj', espressopp.analysis.PotentialEnergy(system, vl_interaction))
    system_analysis.add_observable('lj-tab', espressopp.analysis.PotentialEnergy(system, cg_vl_interaction))
    for lb, interHarmonic in bondedinteractions.iteritems():
        system_analysis.add_observable('bond_%d' % lb, espressopp.analysis.PotentialEnergy(system, interHarmonic), False)
    for lb, interAngularHarmonic in angleinteractions.iteritems():
        system_analysis.add_observable('angle_%d' % lb, espressopp.analysis.PotentialEnergy(system, interAngularHarmonic), False)
    for lb, interDihHarmonic in dihedralinteractions.iteritems():
       system_analysis.add_observable('dihedral_%d' % lb, espressopp.analysis.PotentialEnergy(system, interDihHarmonic), False)
    for lb, interaction14 in pairinteractions.iteritems():
        system_analysis.add_observable('lj-14_%d' % lb, espressopp.analysis.PotentialEnergy(system, interaction14), False)
    system_analysis.add_observable('P', pressure_comp, True)
    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, args.energy_collect)
    integrator.addExtension(ext_analysis)
    print('Configured system analysis')

    print('Save trajectory to {}{}'.format(args.output_prefix, filename))
    ps = espressopp.analysis.PyStore(
        system, '{}{}'.format(args.output_prefix, filename), group_name=h5md_group,
        static_box=False, author='Jakub Krajniak', email='jkrajniak@gmail.com')

    print('Reset total velocity')
    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    ps.dump(0, 0)

    k_trj_collect = args.trj_collect / integrator_step
    k_energy_collect = args.energy_collect / integrator_step

    print('Running simulation for {} steps'.format(sim_step*integrator_step))
    print('Collect trajectory every {} step'.format(k_trj_collect))
    print('Collect energy every {} step'.format(k_energy_collect))
    system_analysis.info()

    for k in range(sim_step):
        time0_ = time.time()
        integrator.run(integrator_step)
        if k_energy_collect > 0 and k % k_energy_collect == 0:
            system_analysis.info()
        if k_trj_collect > 0 and k % k_trj_collect == 0:
            ps.dump(k*integrator_step, k*integrator_step*dt)
        if k_trj_collect > 0 and k % 100 == 0:
            ps.flush()

    ps.dump(sim_step*integrator_step, sim_step*integrator_step*dt)
    ps.close()

    print('finished!')
    print('total time: {}'.format(time.time()-time0))
    espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())


if __name__ == '__main__':
    main()
