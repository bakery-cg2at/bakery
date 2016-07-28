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

import collections
import os

import espressopp  # noqa

__doc__ = 'The tools for the simulation.'


def setSystemAnalysis(system, integrator, args, interval, filename_suffix=None,
                      dynamic_res=None, particle_groups=None):
    """Sets system analysis routine"""
    if filename_suffix is None:
        filename_suffix = ''
    try:
        alpha = args.alpha
        energy_file = '{}_energy_{}_{}{}.csv'.format(args.output_prefix, alpha, args.rng_seed, filename_suffix)
    except AttributeError:
        energy_file = '{}_energy_{}{}.csv'.format(args.output_prefix, args.rng_seed, filename_suffix)
    print('Energy saved to: {}'.format(energy_file))
    system_analysis = espressopp.analysis.SystemMonitor(
        system,
        integrator,
        espressopp.analysis.SystemMonitorOutputCSV(energy_file))
    temp_comp = espressopp.analysis.Temperature(system)
    system_analysis.add_observable('T', temp_comp)

    if particle_groups is not None:
        for label, pg in particle_groups.items():
            system_analysis.add_observable(
                'T-{}'.format(label),
                espressopp.analysis.TemperatureOnGroup(system, pg))

    system_analysis.add_observable('Ekin', espressopp.analysis.KineticEnergy(system, temp_comp))
    if dynamic_res is not None:
        system_analysis.add_observable(
            'res', espressopp.analysis.Resolution(system, dynamic_res))

    try:
        system_info_filter = args.system_info_filter.split(',')
    except AttributeError:
        system_info_filter = None

    for label, interaction in sorted(system.getAllInteractions().items()):
        show_in_system_info = True
        if system_info_filter:
            show_in_system_info = False
            for v in system_info_filter:
                if v in label:
                    show_in_system_info = True
                    break
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction),
            show_in_system_info)

    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, interval)
    integrator.addExtension(ext_analysis)
    return ext_analysis, system_analysis


def setLennardJonesInteractions(system, input_conf, verletlist, cutoff, nonbonded_params=None,  # NOQA
                                hadress=False, ftpl=None, interaction=None, table_groups=[]):
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
            if (pi.get('atnum') in table_groups and pj.get('atnum') in table_groups) or (
                pi.get('atname') in table_groups and pj.get('atname') in table_groups):
                type_pairs.add(tuple(sorted([type_1, type_2])))
            elif pi['particletype'] != 'V' and pj['particletype'] != 'V':
                type_pairs.add(tuple(sorted([type_1, type_2])))
            else:
                print('Skip {}-{} for LJ potential'.format(type_1, type_2))
    type_pairs = sorted(type_pairs)

    if not type_pairs:
        return None

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
            print "Setting LJ interaction for", type_1, type_2, "to sig ", sig, "eps", eps, "cutoff", cutoff
            ljpot = espressopp.interaction.LennardJones(epsilon=eps, sigma=sig, shift='auto',
                                                        cutoff=cutoff)
            if ftpl:
                interaction.setPotentialAT(type1=type_1, type2=type_2, potential=ljpot)
            else:
                interaction.setPotential(type1=type_1, type2=type_2, potential=ljpot)
    return interaction


def setTabulatedInteractions(system, atomtypeparams, vl, cutoff, interaction=None, table_groups=[]):
    """Sets tabulated potential for types that has particletype set to 'V'."""
    spline_type = 2
    if interaction is None:
        interaction = espressopp.interaction.VerletListTabulated(vl)

    type_pairs = set()
    for type_1, v1 in atomtypeparams.iteritems():
        for type_2, v2 in atomtypeparams.iteritems():
            if v1.get('particletype', 'A') == 'V' and v2.get('particletype', 'A') == 'V':
                type_pairs.add(tuple(sorted([type_1, type_2])))
            elif (v1.get('atnum') in table_groups and v2.get('atnum') in table_groups) or (
                v1.get('atname') in table_groups and v2.get('atname') in table_groups):
                type_pairs.add(tuple(sorted([type_1, type_2])))
            else:
                print('Skip {}-{} for tabulated potential'.format(type_1, type_2))
    if not type_pairs:
        return None
    for type_1, type_2 in type_pairs:
        print('Set tabulated potential {}-{}'.format(type_1, type_2))
        name_1 = atomtypeparams[type_1]['atnum']
        name_2 = atomtypeparams[type_2]['atnum']
        table_name = '{}-{}.espp.pot'.format(name_1, name_2)
        orig_table_name = 'table_{}_{}.xvg'.format(name_1, name_2)
        if not os.path.exists(table_name):
            espressopp.tools.convert.gromacs.convertTable(orig_table_name, table_name)
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
                tmp.append(0.0)
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


def setBondedInteractions(system, input_conf, force_static=False, only_at=False, only_cg=None):
    ret_list = {}
    bonds = input_conf.bondtypes
    bondtypeparams = input_conf.bondtypeparams

    for (bid, cross_bonds), bondlist in bonds.iteritems():
        b1 = bondlist[0][0]
        is_cg = input_conf.atomtypeparams[input_conf.types[b1-1]]['particletype'] == 'V'

        if only_at and is_cg:
            continue

        if only_cg is not None and only_cg and not is_cg:
            continue

        fpl = espressopp.FixedPairList(system.storage)

        fpl.addBonds(bondlist)
        if not cross_bonds or force_static:
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
        if params['sig'] > 0.0 and params['eps'] > 0.0:
            is_cg = input_conf.atomtypeparams[
                        input_conf.types[pair_list[0][0] - 1]]['particletype'] == 'V'
            if is_cg or ftpl is None:
                fpl = espressopp.FixedPairList(system.storage)
            else:
                fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
            fpl.addBonds(pair_list)

            if not cross_bonds:
                is_cg = None

            print ('Pair interaction', params, ' num pairs:', len(pair_list),
                   'sig=', params['sig'], params['eps'])
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
            system.addInteraction(interaction, 'lj-14_{}{}'.format(pid, '_cross' if cross_bonds else ''))


def setAngleInteractions(system, input_conf, force_static=False, only_at=False, only_cg=None):
    ret_list = {}
    angletypeparams = input_conf.angletypeparams
    angles = input_conf.angletypes

    for (aid, cross_angles), anglelist in angles.iteritems():
        b1 = anglelist[0][0]
        is_cg = input_conf.atomtypeparams[input_conf.types[b1-1]]['particletype'] == 'V'

        if only_at and is_cg:
            continue

        if only_cg is not None and only_cg and not is_cg:
            continue

        fpl = espressopp.FixedTripleList(system.storage)

        fpl.addTriples(anglelist)
        if not cross_angles or force_static:
            is_cg = None
        angleinteraction = angletypeparams[aid].createEspressoInteraction(system, fpl, is_cg=is_cg)
        if angleinteraction:
            system.addInteraction(angleinteraction, 'angle_{}{}'.format(
                aid, '_cross' if cross_angles else ''))
            ret_list.update({(aid, cross_angles): angleinteraction})
    return ret_list


def setDihedralInteractions(system, input_conf, force_static=False, only_at=False,  only_cg=None):
    ret_list = {}
    dihedrals = input_conf.dihedraltypes
    dihedraltypeparams = input_conf.dihedraltypeparams

    for (did, cross_dih), dihedrallist in dihedrals.iteritems():
        b1 = dihedrallist[0][0]
        is_cg = input_conf.atomtypeparams[input_conf.types[b1-1]]['particletype'] == 'V'

        if only_at and is_cg:
            continue

        if only_cg is not None and only_cg and not is_cg:
            continue

        fpl = espressopp.FixedQuadrupleList(system.storage)
        fpl.addQuadruples(dihedrallist)
        if not cross_dih or force_static:
            is_cg = None
        dihedralinteraction = dihedraltypeparams[did].createEspressoInteraction(
            system, fpl, is_cg=is_cg)
        if dihedralinteraction:
            system.addInteraction(dihedralinteraction, 'dihedral_{}{}'.format(
                did, '_cross' if cross_dih else ''))
            ret_list.update({(did, cross_dih): dihedralinteraction})
    return ret_list


def saveInteractions(system, output_filename):
    number_of_interactions = system.getNumberOfInteractions()
    print('Number of interactions: {}'.format(number_of_interactions))

    data = {}
    for interaction_name, interaction in system.getAllInteractions().items():
        print('Interaction {}'.format(interaction_name))
        tmp_data = {'class_str': str(interaction)}
        tmp_data['params'] = interaction.getParams()
        fixed_list = interaction.getFixedList()
        if fixed_list is not None:
            tmp_data['fixed_list'] = fixed_list.getList()
        data[interaction_name] = tmp_data

    with open(output_filename, 'wb') as output_file:
        import cPickle
        cPickle.dump(data, output_file)
        print('Saved in {}'.format(output_filename))