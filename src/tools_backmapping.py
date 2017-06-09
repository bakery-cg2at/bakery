#! /usr/bin/env python
#
# Copyright (c) 2016,2017 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import espressopp
import tools_sim as tools
import gromacs_topology
import cPickle


def setupSinglePhase(system, args, input_conf, at_particle_ids, cg_particle_ids, table_groups=[]):
    exclusionlistAT = [p for p in input_conf.exclusions
                       if p[0] in at_particle_ids and p[1] in at_particle_ids]
    exclusionlistCG = [p for p in input_conf.exclusions
                       if p[0] in cg_particle_ids and p[1] in cg_particle_ids]
    print('Excluded pairs for LJ interaction (AT): {}'.format(len(exclusionlistAT)))
    print('Excluded pairs for LJ interaction (CG): {}'.format(len(exclusionlistCG)))
    verletlistAT = espressopp.VerletListHybridAT(
        system, cutoff=args.lj_cutoff, exclusionlist=exclusionlistAT)

    verletlistCG = espressopp.VerletListHybridCG(
        system, cutoff=args.cg_cutoff, exclusionlist=exclusionlistCG)

    lj_interaction = espressopp.interaction.VerletListHybridLennardJones(verletlistAT, False)
    if args.cap_force_lj:
        print('Defined max_force for LJ potential'.format(args.cap_force_lj))
        lj_interaction.max_force = args.cap_force_lj
    lj_interaction = tools.setLennardJonesInteractions(
        system, input_conf, verletlistAT, args.lj_cutoff,
        input_conf.nonbond_params,
        interaction=lj_interaction,
        table_groups=table_groups)

    coulomb_interaction = espressopp.interaction.VerletListHybridReactionFieldGeneralized(
        verletlistAT, False)
    #if args.cap_force_lj:
    #    coulomb_interaction.max_force = args.cap_force_lj
    if args.coulomb_cutoff > 0.0:
        coulomb_interaction = gromacs_topology.setCoulombInteractions(
            system, verletlistAT, args.coulomb_cutoff, input_conf.atomtypeparams,
            epsilon1=args.coulomb_epsilon1,
            epsilon2=args.coulomb_epsilon2, kappa=args.coulomb_kappa,
            interaction=coulomb_interaction)
    else:
        coulomb_interaction = None
    ret_list = tools.setBondedInteractions(
        system, input_conf)
    tools.setAngleInteractions(
        system, input_conf)
    tools.setDihedralInteractions(
        system, input_conf)
    # tools.setPairInteractions(
    #     system, input_conf, args.lj_cutoff, args.coulomb_cutoff)
    tab_cg_interaction = espressopp.interaction.VerletListHybridTabulated(verletlistCG, True)
    #if args.cap_force_lj:
    #    tab_cg_interaction.max_force = args.cap_force_lj
    tab_cg = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlistCG,
        cutoff=args.cg_cutoff,
        interaction=tab_cg_interaction,
        table_groups=table_groups)
    if lj_interaction is not None:
        system.addInteraction(lj_interaction, 'lj')
    if coulomb_interaction is not None:
        system.addInteraction(coulomb_interaction, 'coulomb')
    if tab_cg is not None:
        system.addInteraction(tab_cg, 'tab-cg')

    # Save interactions
    if args.save_interactions:
        tools.saveInteractions(system, '{}_{}_interactions.pck'.format(args.output_prefix, args.rng_seed))

    return verletlistAT, verletlistCG


def setupFirstPhase(system, args, input_conf, at_particle_ids, cg_particle_ids):
    """Two phase. First keep cg interactions static and scale only bonded."""
    # First remove all interactions and define again.
    print('Setup first phase in two-phase backmapping scheme')
    number_of_interactions = system.getNumberOfInteractions()-1
    for interaction_id in range(number_of_interactions, -1, -1):
        system.removeInteraction(interaction_id)

    # Reset lambda for all particles.
    for pid in at_particle_ids:
        system.storage.modifyParticle(pid, 'lambda_adr', 0.0)
    for pid in cg_particle_ids:
        system.storage.modifyParticle(pid, 'lambda_adr', 0.0)

    system.storage.decompose()

    exclusionlistCG = [p for p in input_conf.exclusions
                       if p[0] in cg_particle_ids and p[1] in cg_particle_ids]
    print('Excluded pairs for LJ interaction (CG): {}'.format(len(exclusionlistCG)))

    verletlistCG = espressopp.VerletListHybridCG(
        system,
        cutoff=args.cg_cutoff,
        exclusionlist=exclusionlistCG)

    tools.setBondedInteractions(
        system, input_conf)
    tools.setAngleInteractions(
        system, input_conf)
    tools.setDihedralInteractions(
        system, input_conf)
    tab_cg = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlistCG,
        cutoff=args.cg_cutoff,
        interaction=espressopp.interaction.VerletListTabulated(verletlistCG))
    if tab_cg is not None:
        system.addInteraction(tab_cg, 'tab-cg')

    if args.save_interactions:
        tools.saveInteractions(system, '{}_{}_phase_one_interactions.pck'.format(args.output_prefix, args.rng_seed))

    return verletlistCG


def setupSecondPhase(system, args, input_conf, at_particle_ids, cg_particle_ids):
    """Two phase. First keep cg interactions static and scale only bonded."""
    print('Setup second phase in two-phase backmapping scheme')
    # First remove all interactions and define again.
    number_of_interactions = system.getNumberOfInteractions()-1
    for interaction_id in range(number_of_interactions, -1, -1):
        system.removeInteraction(interaction_id)

    # Reset lambda for all particles.
    for pid in at_particle_ids:
        system.storage.modifyParticle(pid, 'lambda_adr', 0.0)
    for pid in cg_particle_ids:
        system.storage.modifyParticle(pid, 'lambda_adr', 0.0)
    system.storage.decompose()

    # Reset all interactions.
    exclusionlistAT = [p for p in input_conf.exclusions
                       if p[0] in at_particle_ids and p[1] in at_particle_ids]
    exclusionlistCG = [p for p in input_conf.exclusions
                       if p[0] in cg_particle_ids and p[1] in cg_particle_ids]
    print('Number of CG particles: {}'.format(len(cg_particle_ids)))
    print('Number of AT particles: {}'.format(len(at_particle_ids)))
    print('Excluded pairs for LJ interaction (AT): {}'.format(len(exclusionlistAT)))
    print('Excluded pairs for LJ interaction (CG): {}'.format(len(exclusionlistCG)))
    verletlistAT = espressopp.VerletListHybridAT(
        system, cutoff=args.lj_cutoff, exclusionlist=exclusionlistAT)

    verletlistCG = espressopp.VerletListHybridCG(
        system, cutoff=args.cg_cutoff, exclusionlist=exclusionlistCG)

    lj_interaction = espressopp.interaction.VerletListHybridLennardJones(
        verletlistAT, False)
    lj_interaction = tools.setLennardJonesInteractions(
        system, input_conf, verletlistAT, args.lj_cutoff,
        input_conf.nonbond_params,
        interaction=lj_interaction)
    if args.coulomb_cutoff > 0.0:
        coulomb_interaction = gromacs_topology.setCoulombInteractions(
            system, verletlistAT, args.coulomb_cutoff, input_conf.atomtypeparams,
            epsilon1=args.coulomb_epsilon1,
            epsilon2=args.coulomb_epsilon2, kappa=args.coulomb_kappa,
            interaction=espressopp.interaction.VerletListHybridReactionFieldGeneralized(
                verletlistAT, False))
    else:
        coulomb_interaction = None
    # tools.setPairInteractions(
    #    system, input_conf, args.lj_cutoff, args.coulomb_cutoff)
    tab_cg = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlistCG,
        cutoff=args.cg_cutoff,
        interaction=espressopp.interaction.VerletListHybridTabulated(
            verletlistCG, True))
    if lj_interaction is not None:
        system.addInteraction(lj_interaction, 'lj')
    if coulomb_interaction is not None:
        system.addInteraction(coulomb_interaction, 'coulomb')
    if tab_cg is not None:
        system.addInteraction(tab_cg, 'lj-tab')

    # The bonded terms are already correct, not scale with lambda.
    ret_list = tools.setBondedInteractions(
        system, input_conf, force_static=True, only_at=True)

    tools.setAngleInteractions(
        system, input_conf, force_static=True, only_at=True)
    tools.setDihedralInteractions(
        system, input_conf, force_static=True, only_at=True)

    if args.save_interactions:
        tools.saveInteractions(system, '{}_{}_phase_two_interactions.pck'.format(args.output_prefix, args.rng_seed))

    return verletlistAT, verletlistCG
