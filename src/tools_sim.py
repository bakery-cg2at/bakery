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
import ast
import collections
import math
import os
import sys

import espressopp  # noqa
import numpy

import random
import tools

__doc__ = 'The tools for the simulation.'


Molecule = collections.namedtuple('Molecule', ['pid', 'pos', 'mass', 'type'])


def replicate_list(N_molecules, N_single, molecule_list, shift=0):
    return [
        map(lambda x: shift+x+(n*N_single), z)
        for n in range(N_molecules) for z in molecule_list
        ]


def pdbread(filename, scale_factor=0.1):
    """Reads the pdb file.

    Warning: currently it is a very basic implementation. Reads only particle
    position and if exists the CONECT section.

    The position and box size is expressed in the Anstrom units.

    Args:
    filename: The input pdb file.
    scale_factor: The multiplicator of the numerical value.

    Returns:
    The optional configuration of box.
    The list of lists with particles grouped by the residue sequence number.
      Each particle is defined by the tuple with properties 'pid' and 'pos'
    The optional list of bonds.
    """
    atoms = []
    atom_bonds = collections.defaultdict(set)
    box = None
    pdb_file = open(filename, 'r')

    file_content = pdb_file.readlines()

    # Following http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
    for line in file_content:
        if line.startswith('CRYST1'):  # crystal
            box = espressopp.Real3D(
                float(line[6:15])*scale_factor,
                float(line[15:24])*scale_factor,
                float(line[24:33])*scale_factor
                )
        elif line.startswith('ATOM') or line.startswith('HETATM'):
            atom_id = int(line[6:11].strip())
            atom_name = line[12:16].strip()
            chain_idx = line[22:26].strip()
            pos = espressopp.Real3D(
                float(line[30:38])*scale_factor,
                float(line[38:46])*scale_factor,
                float(line[46:54])*scale_factor)
            atoms.append((atom_id, chain_idx, atom_name, pos))
        elif line.startswith('CONECT'):
            # http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_69.html
            bonded_atoms = map(int, line[6:].split())
            if bonded_atoms:
                at0 = bonded_atoms[0]
                at1 = set(bonded_atoms[1:])
                atom_bonds[at0].update(at1)

    bonds = set()
    for at0, ats in atom_bonds.iteritems():
        for at in ats:
            bond = (at0, at)
            if not ((at0, at) in bonds or (at, at0) in bonds):
                bonds.add(bond)

    return box, atoms, list(bonds)


def groread(filename, scale_factor=1.0):
    """Reads the .gro file and return the atom list.

    Returns:
        The dict with atoms (key: atom_id, value: atom object).
    """

    atoms = []
    input_file = open(filename, 'r')
    content = input_file.readlines()

    number_of_atoms = int(content[1])

    for line in content[2:number_of_atoms + 2]:
        chain_idx = int(line[0:5].strip())
        # chain_name = line[5:10].strip()
        at_name = line[10:15].strip()
        at_id = int(line[15:20].strip())
        # Nedd to rescale.
        pos_x = float(line[20:28].strip()) * scale_factor
        pos_y = float(line[28:36].strip()) * scale_factor
        pos_z = float(line[36:44].strip()) * scale_factor
        atoms.append((at_id, chain_idx, at_name, espressopp.Real3D(pos_x, pos_y, pos_z)))

    # Reads the box size, the last line.
    box = espressopp.Real3D(numpy.array(
        map(float, filter(None, content[number_of_atoms + 2].split(' ')))
        ) * scale_factor)

    return box, atoms, []


def growrite(atoms, box, file_name):
    """Writes the content to the output file.

    Args:
        atoms: The dict with atoms.
        box: The tuple with box description.
        file_name: The new file name, otherwise the old one will be used.
    """

    output = ['XXX of molecules']
    # Puts the number of atoms
    output.append('%d' % len(atoms))
    # Puts the definition of the atoms, fixed format.
    fmt = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
    for at_id, chain_idx, at_name, pos in atoms:
        output.append(fmt % (
            chain_idx,
            'XXX',
            at_name,
            at_id,
            pos[0],
            pos[1],
            pos[2]
            ))

    output.append('%f %f %f' % tuple(box))
    output_file = open(file_name, 'w')
    output_file.writelines('\n'.join(output))
    output_file.close()


def readtrj(filename):
    return {
        'gro': groread,
        'pdb': pdbread
    }[filename.split('.')[-1]](filename)


def warmup(system,
           integrator,
           number,
           potential_pairs=None,
           total_force_capping=1000000.0):
    """Warmup the system.

    Args:
        system: The system object.
        integrator: The integrator object.
        number: The number of steps of the warm up phase.
        potential_pairs: The list of tuples with the particle types potentials.
        total_force_capping: The maximum force capping.
    """
    if potential_pairs is None:
        potential_pairs = [(0, 0)]

    print('\nStarting warmup')

    org_dt = integrator.dt
    integrator.dt = pow(10, -4)

    N = 100
    warmup_steps = number / N
    d_force = total_force_capping / warmup_steps

    print('warmup with dt = %e' % integrator.dt)
    print('max_force = %d' % total_force_capping)
    print('number = %d' % number)
    print('warmup time = %e s' % (warmup_steps*integrator.dt))
    print('d_force = %e' % d_force)

    potentials = {}
    interaction = system.getInteractionByLabel('lj')
    for pair in potential_pairs:
        pot = interaction.getPotential(*pair)
        potentials[pair] = [
            pot,
            pot.sigma,
            pot.epsilon,
            pot.epsilon / warmup_steps,
            pot.sigma / warmup_steps
            ]

    force_capping = espressopp.integrator.CapForce(system, d_force)
    integrator.addExtension(force_capping)
    info(system, integrator)

    for k in range(1, warmup_steps+1):
        # Update the sigma and the epsilon.
        for pair, pot in potentials.iteritems():
            pot[0].epsilon = k*pot[3]
            pot[0].sigma = k*pot[4]
            interaction.setPotential(pair[0], pair[1], pot[0])
        # e_pots.append(EPot.compute())
        # temps.append(T.compute())
        # print k, numpy.std(e_pots), numpy.std(temps)
        integrator.run(N)
        info(system, integrator)
        force_capping.setAbsCapForce(k*d_force)

    # Restore parameters.
    for pair, pot in potentials.iteritems():
        pot[0].sigma = pot[1]
        pot[0].epsilon = pot[2]
        interaction.setPotential(pair[0], pair[1], pot[0])

    force_capping.disconnect()

    print('Force capping disconnected')

    for k in range(2*warmup_steps):
        info(system, integrator)
        integrator.run(N)

    integrator.dt = org_dt
    integrator.step = 0
    print("warmup finished")


# Compute the angle distribution
def unit_vector(vector):
    """Returns the unit vector.

    Args:
    vector: The input vector.

    Returns:
    Returns the normalized vector.
    """
    return vector / numpy.linalg.norm(vector)


def info(system, integrator, per_atom=False):
    """Display some information during the simulation.

    Args:
        system: The system object.
        integrator: The integrator object.
        per_atom: Compute per atom.

    Returns:
        The list with step, temperature, pressure, pressure_xy
        kinetic energy,
        potential energy for each of the interaction,
        total potential, total energy.
    """

    NPart = espressopp.analysis.NPart(system).compute()
    T = espressopp.analysis.Temperature(system).compute() / system.kb
    P = espressopp.analysis.Pressure(system).compute()
    Pij = espressopp.analysis.PressureTensor(system).compute()
    step = integrator.step
    Ek = (3.0/2.0) * NPart * T
    Etotal = 0.0

    if any(map(math.isnan, [T, P, Pij[3]])) or any(map(math.isinf, [T, P, Pij[3]])):
        raise ValueError('Temperature, pressure or pressure tensor is not valid')

    if per_atom:
        data = [step, step*integrator.dt, T, P, Pij[3], Ek/NPart]
        tot = '%5d %10.4f %10.6f %10.6f %10.6f %12.8f' % tuple(data)
    else:
        data = [step, step*integrator.dt, T, P, Pij[3], Ek]
        tot = '%5d %10.4f %10.6f %10.6f %10.6f %12.3f' % tuple(data)
    header = ''
    for name, k in system.getInterationLabels().iteritems():
        e = system.getInteraction(k).computeEnergy()
        if math.isnan(e) or math.isinf(e):
            raise ValueError('The value for {} is not valid.'.format(name))
        Etotal += e
        if per_atom:
            tot += ' %12.8f' % (e/NPart)
            data.append(e/NPart)
            header += '     %s%i/N    ' % (name, k)
        else:
            tot += ' %12.3f' % e
            data.append(e)
            header += '      %s%i     ' % (name, k)

    if per_atom:
        tot += ' %12.8f' % (Etotal/NPart)
        tot += ' %12.8f' % (Etotal/NPart + Ek/NPart)
        data.append(Etotal/NPart)
        data.append(Etotal/NPart + Ek/NPart)
        header += '   epot/N  ' + '   etotal/N  '
    else:
        tot += ' %12.8f' % (Etotal)
        tot += ' %12.3f' % (Etotal + Ek)
        data.append(Etotal)
        data.append(Etotal + Ek)
        header += '   epot/N  ' + '   etotal/N  '

    tot += ' %12.8f\n' % system.bc.boxL[0]
    header += '    boxL     \n'
    if step == 0:
        if per_atom:
            sys.stdout.write(' step      dt     T          P        Pxy         ekin/N  ' + header)
        else:
            sys.stdout.write(' step      dt     T          P        Pxy          ekin   ' + header)
        sys.stdout.write(tot)

    return data


def final_info(system, integrator, vl, start_time, end_time):
    NPart = espressopp.analysis.NPart(system).compute()
    espressopp.tools.timers.show(integrator.getTimers(), precision=3)
    sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
    sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(NPart)))
    sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
    sys.stdout.write('Integration steps = %d\n' % integrator.step)
    sys.stdout.write('CPUs = %i CPU time per CPU = %.5f\n' % (
        espressopp.MPI.COMM_WORLD.size, end_time - start_time))


def warmup_capped(system, integrator, verletList, rc, potential_matrix,
                  sigma_start, warmup_loops, warmup_steps):
    LJ = espressopp.interaction.VerletListLennardJones(verletList)
    for type1, type2, sigma, epsilon in potential_matrix:
        LJ.setPotential(
            type1=type1,
            type2=type2,
            potential=espressopp.interaction.LennardJones(
                sigma=sigma,
                epsilon=epsilon,
                cutoff=rc
            ))
    system.addInteraction(LJ, 'lj_warmup')

    old_dt = integrator.dt
    integrator.dt = 0.0001

    espressopp.tools.analyse.info(system, integrator, kb=system.kb)

    print('Running with capped potential, increasing sigma, loops {}, steps {}'.format(
        warmup_loops, warmup_steps))
    # Run system with capped potentials, thermostat and increasing LJ epsilon
    for k in range(warmup_loops):
        for type1, type2, sigma, epsilon in potential_matrix:
            LJ.setPotential(
                type1=type1,
                type2=type2,
                potential=espressopp.interaction.LennardJones(
                    sigma=sigma_start + (sigma-sigma_start)*k*1.0/(warmup_loops-1),
                    epsilon=epsilon,
                    cutoff=rc
                ))
        integrator.run(warmup_steps)
        if integrator.step % 100 == 0:
            espressopp.tools.analyse.info(system, integrator, kb=system.kb)

    # Remove LJ Capped potential
    system.removeInteractionByLabel('lj_warmup')
    integrator.dt = old_dt


def setLennardJonesInteractions(system, input_conf, verletlist, cutoff, nonbonded_params=None,  # NOQA
                                hadress=False, ftpl=None, interaction=None):
    """ Set lennard jones interactions which were read from gromacs based on the atomypes"""
    defaults = input_conf.defaults
    atomtypeparams = input_conf.atomtypeparams
    if interaction is None:
        if ftpl:
            if hadress:
                interaction = espressopp.interaction.VerletListHadressLennardJones(verletlist, ftpl)
            else:
                interaction = espressopp.interaction.VerletListAdressLennardJones(verletlist, ftpl)
        else:
            interaction = espressopp.interaction.VerletListLennardJones(verletlist)

    if nonbonded_params is None:
        nonbonded_params = {}

    combinationrule = int(defaults['combinationrule'])
    print "Setting up Lennard-Jones interactions"

    type_pairs = set()
    for type_1, pi in atomtypeparams.iteritems():
        for type_2, pj in atomtypeparams.iteritems():
            if pi['particletype'] != 'V' and pj['particletype'] != 'V':
                type_pairs.add(tuple(sorted([type_1, type_2])))
    type_pairs = sorted(type_pairs)

    print('Number of pairs: {}'.format(len(type_pairs)))
    for type_1, type_2 in type_pairs:
        pi = atomtypeparams[type_1]
        pj = atomtypeparams[type_2]
        if pi['particletype'] == 'V' or pj['particletype'] == 'V':
            print('Skip {}-{}'.format(type_1, type_2))
            continue
        param = nonbonded_params.get((type_1, type_2))
        if param:
            print 'Using defined non-bonded cross params', param
            sig, eps = param['sig'], param['eps']
        else:
            sig_1, eps_1 = float(pi['sig']), float(pi['eps'])
            sig_2, eps_2 = float(pj['sig']), float(pj['eps'])
            if combinationrule == 2:
                sig = 0.5*(sig_1 + sig_2)
                eps = (eps_1*eps_2)**(1.0/2.0)
            else:
                sig = (sig_1*sig_2)**(1.0/2.0)
                eps = (eps_1*eps_2)**(1.0/2.0)
        if sig > 0.0 and eps > 0.0:
            print ("Setting LJ interaction for", type_1, type_2, "to sig ", sig, "eps",
                   eps, "cutoff", cutoff)
            ljpot = espressopp.interaction.LennardJones(epsilon=eps, sigma=sig, shift='auto',
                                                        cutoff=cutoff)
            if ftpl:
                interaction.setPotentialAT(type1=type_1, type2=type_2, potential=ljpot)
            else:
                interaction.setPotential(type1=type_1, type2=type_2, potential=ljpot)
    return interaction


def setTabulatedInteractions(system, atomtypeparams, vl, cutoff, interaction=None, ftpl=None):
    """Sets tabulated potential for types that has particletype set to 'V'."""
    spline_type = 2
    if interaction is None:
        if ftpl:
            interaction = espressopp.interaction.VerletListAdressTabulated(vl, ftpl)
        else:
            interaction = espressopp.interaction.VerletListTabulated(vl)
    type_pairs = set()
    for type_1, v1 in atomtypeparams.iteritems():
        for type_2, v2 in atomtypeparams.iteritems():
            if v1.get('particletype', 'A') == 'V' and v2.get('particletype', 'A') == 'V':
                type_pairs.add(tuple(sorted([type_1, type_2])))
    for type_1, type_2 in type_pairs:
        print('Set tabulated potential {}-{}'.format(type_1, type_2))
        name_1 = atomtypeparams[type_1]['atnum']
        name_2 = atomtypeparams[type_2]['atnum']
        table_name = '{}-{}.espp.pot'.format(name_1, name_2)
        orig_table_name = 'table_{}_{}.xvg'.format(name_1, name_2)
        if not os.path.exists(table_name):
            espressopp.tools.convert.gromacs.convertTable(orig_table_name, table_name)
        if ftpl:
            interaction.setPotentialCG(
                type1=type_1,
                type2=type_2,
                potential=espressopp.interaction.Tabulated(
                    itype=spline_type,
                    filename=table_name,
                    cutoff=cutoff))
        else:
            interaction.setPotential(
                type1=type_1,
                type2=type_2,
                potential=espressopp.interaction.Tabulated(
                    itype=spline_type,
                    filename=table_name,
                    cutoff=cutoff))
    return interaction


def genParticleList(input_conf, use_velocity=False, use_charge=False, adress=False):  #NOQA
    """Generates particle list
    Args:
        input_conf: The tuple generate by read method.
        use_velocity: If set to true then velocity will be read.
        use_charge: If set to true then charge will be read.
        adress: If set to true then adress_tuple will be generated.
    Returns:
        List of property names and particle list.
    """
    props = ['id', 'type', 'pos', 'res_id']
    use_mass = bool(input_conf.masses)
    use_velocity = use_velocity and bool(input_conf.vx)
    use_charge = use_charge and bool(input_conf.charges)
    if use_mass:
        props.append('mass')
    if use_velocity:
        props.append('v')
    if use_charge:
        props.append('q')

    Particle = collections.namedtuple('Particle', props)
    particle_list = []
    num_particles = len(input_conf.types)
    if adress:
        props.append('adrat')   # Set to 1 if AT particle otherwise 0
        props.append('lambda_adr')
        props.append('vp')
        Particle = collections.namedtuple('Particle', props)
        adress_tuple = []
        tmptuple = []
        for pid in range(num_particles):
            atom_type = input_conf.types[pid]
            particle_type = input_conf.atomtypeparams[atom_type]['particletype']
            tmp = [pid+1,
                   atom_type,
                   espressopp.Real3D(input_conf.x[pid], input_conf.y[pid], input_conf.z[pid]),
                   input_conf.res_ids[pid]]
            if use_mass:
                tmp.append(input_conf.masses[pid])
            if use_velocity:
                tmp.append(espressopp.Real3D(
                    input_conf.vx[pid],
                    input_conf.vy[pid],
                    input_conf.vz[pid]))
            if use_charge:
                tmp.append(input_conf.charges[pid])
            if particle_type == 'V':
                tmp.append(0)  # adrat
                tmp.append(1.0)
                tmp.append(True)
                if tmptuple != []:
                    adress_tuple.append(tmptuple[:])
                tmptuple = [pid+1]
            else:
                tmp.append(1)  # adrat
                tmp.append(0.0)
                tmp.append(False)
                tmptuple.append(pid+1)
            particle_list.append(Particle(*tmp))
        # Set Adress tuples
        adress_tuple.append(tmptuple[:])
        return props, particle_list, adress_tuple
    else:
        for pid in range(num_particles):
            tmp = [pid+1,
                   input_conf.types[pid],
                   espressopp.Real3D(input_conf.x[pid], input_conf.y[pid], input_conf.z[pid]),
                   input_conf.res_ids[pid]]
            if use_mass:
                tmp.append(input_conf.masses[pid])
            if use_velocity:
                tmp.append(espressopp.Real3D(
                    input_conf.vx[pid],
                    input_conf.vy[pid],
                    input_conf.vz[pid]))
            if use_charge:
                tmp.append(input_conf.charges[pid])
            particle_list.append(Particle(*tmp))
        return props, particle_list


def setBondedInteractions(system, input_conf, ftpl=None):
    ret_list = {}
    bonds = input_conf.bondtypes
    bondtypeparams = input_conf.bondtypeparams

    for (bid, cross_bonds), bondlist in bonds.iteritems():
        b1 = bondlist[0][0]
        is_cg = input_conf.atomtypeparams[input_conf.types[b1-1]]['particletype'] == 'V'

        fpl = espressopp.FixedPairList(system.storage)

        fpl.addBonds(bondlist)
        if not cross_bonds:
            is_cg = None
        bdinteraction = bondtypeparams[bid].createEspressoInteraction(system, fpl, is_cg=is_cg)
        if bdinteraction:
            system.addInteraction(bdinteraction, 'bond_{}{}'.format(
                bid, '_cross' if cross_bonds else ''))
            ret_list.update({(bid, cross_bonds): bdinteraction})

    return ret_list


def setPairInteractions(system, input_conf, cutoff, ftpl=None):
    ret_list = {}
    pairs = input_conf.pairtypes
    pairtypeparams = input_conf.pairtypeparams
    for (pid, cross_bonds), pair_list in pairs.iteritems():
        params = pairtypeparams[pid]
        is_cg = input_conf.atomtypeparams[
            input_conf.types[pair_list[0][0]-1]]['particletype'] == 'V'
        if is_cg or ftpl is None:
            fpl = espressopp.FixedPairList(system.storage)
        else:
            fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
        fpl.addBonds(pair_list)
        print ('Pair interaction', params, ' num pairs:', len(pair_list),
               'sig=', params['sig'], params['eps'])

        if not cross_bonds:
            is_cg = None

        pot = espressopp.interaction.LennardJones(
            sigma=params['sig'],
            epsilon=params['eps'],
            shift='auto',
            cutoff=cutoff)
        if is_cg is None:
            interaction = espressopp.interaction.FixedPairListLennardJones(system, fpl, pot)
        else:
            interaction = espressopp.interaction.FixedPairListAdressLennardJones(
                system, fpl, pot, is_cg)
        system.addInteraction(interaction, 'lj-14_{}{}'.format(
            pid, '_cross' if cross_bonds else ''))
        ret_list[(pid, cross_bonds)] = interaction
    return ret_list


def setAngleInteractions(system, input_conf, ftpl=None):
    ret_list = {}
    angletypeparams = input_conf.angletypeparams
    angles = input_conf.angletypes

    for (aid, cross_angles), anglelist in angles.iteritems():
        b1 = anglelist[0][0]
        is_cg = input_conf.atomtypeparams[input_conf.types[b1-1]]['particletype'] == 'V'

        if is_cg or ftpl is None:
            fpl = espressopp.FixedTripleList(system.storage)
        else:
            fpl = espressopp.FixedTripleListAdress(system.storage, ftpl)
        fpl.addTriples(anglelist)
        if not cross_angles:
            is_cg = None
        angleinteraction = angletypeparams[aid].createEspressoInteraction(system, fpl, is_cg=is_cg)
        if angleinteraction:
            system.addInteraction(angleinteraction, 'angle_{}{}'.format(
                aid, '_cross' if cross_angles else ''))
            ret_list.update({(aid, cross_angles): angleinteraction})
    return ret_list


def setDihedralInteractions(system, input_conf, ftpl=None):
    ret_list = {}
    dihedrals = input_conf.dihedraltypes
    dihedraltypeparams = input_conf.dihedraltypeparams

    for (did, cross_dih), dihedrallist in dihedrals.iteritems():
        b1 = dihedrallist[0][0]
        is_cg = input_conf.atomtypeparams[input_conf.types[b1-1]]['particletype'] == 'V'

        if is_cg or ftpl is None:
            fpl = espressopp.FixedQuadrupleList(system.storage)
        else:
            fpl = espressopp.FixedQuadrupleListAdress(system.storage, ftpl)
        fpl.addQuadruples(dihedrallist)
        if not cross_dih:
            is_cg = None
        dihedralinteraction = dihedraltypeparams[did].createEspressoInteraction(
            system, fpl, is_cg=is_cg)
        if dihedralinteraction:
            system.addInteraction(dihedralinteraction, 'dihedral_{}{}'.format(
                did, '_cross' if cross_dih else ''))
            ret_list.update({(did, cross_dih): dihedralinteraction})
    return ret_list


def _args_md():
    parser = tools.MyArgParser(description='Runs classical MD simulation',
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
    parser.add_argument('--coulomb_epsilon2', default=80.0, type=float,
                        help='Epsilon_2 for coulomb interactions')
    parser.add_argument('--coulomb_kappa', default=1.0, type=float,
                        help='Kappa paramter for coulomb interactions')
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
