"""The tools for the simulation."""


import collections
import math
import sys

import espressopp  # noqa
import numpy


Molecule = collections.namedtuple('Molecule', ['pid', 'pos', 'mass', 'type'])


class System(espressopp.System):
    """Helper class to manage interactions."""

    def __init__(self, kb=1.0):
        super(System, self).__init__()
        self._interaction_id = 0
        self._interaction_label_id = {}
        self.kb = kb

    def addInteraction(self, interaction, label=None):
        if label is None:
            label = 'e%d' % self._interaction_id
        if label is not None and label in self._interaction_label_id:
            raise ValueError('Interaction with label %s exists', label)
        print('Adding interaction {}: {}'.format(label, interaction))
        super(System, self).addInteraction(interaction)
        if label is not None:
            self._interaction_label_id[label] = self._interaction_id
            self._interaction_id += 1

    def removeInteractionByLabel(self, label):
        interaction_id = self._interaction_label_id[label]
        super(System, self).removeInteraction(interaction_id)
        del self._interaction_label_id[label]
        self._interaction_id -= 1

    def getInteractionByLabel(self, label):
        return super(System, self).getInteraction(self._interaction_label_id[label])

    def getInteractionLabels(self):
        return self._interaction_label_id

    def printInteractions(self):
        label_ids = sorted(self._interaction_label_id.items(), key=lambda x: x[1])
        for label, interaction_id in label_ids:
            print('{} -> e{}'.format(label, interaction_id))


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
