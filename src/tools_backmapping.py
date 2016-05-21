#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import espressopp
import tools_sim as tools
import gromacs_topology


def setupSinglePhase(system, args, input_conf, at_particle_ids, cg_particle_ids):
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

    lj_interaction = espressopp.interaction.VerletListHybridLennardJones(
        verletlistAT, False)
    lj_interaction = tools.setLennardJonesInteractions(
        system, input_conf, verletlistAT, args.lj_cutoff,
        input_conf.nonbond_params,
        interaction=lj_interaction)
    coulomb_interaction = espressopp.interaction.VerletListHybridReactionFieldGeneralized(
        verletlistAT, False)
    coulomb_interaction = gromacs_topology.setCoulombInteractions(
        system, verletlistAT, 0.9, input_conf.atomtypeparams,
        epsilon1=args.coulomb_epsilon1,
        epsilon2=args.coulomb_epsilon2, kappa=args.coulomb_kappa,
        interaction=coulomb_interaction)
    tools.setBondedInteractions(
        system, input_conf)
    tools.setAngleInteractions(
        system, input_conf)
    tools.setDihedralInteractions(
        system, input_conf)
    pair14_interactions = tools.setPairInteractions(
        system, input_conf, args.lj_cutoff)
    tab_cg = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlistCG,
        cutoff=args.cg_cutoff,
        interaction=espressopp.interaction.VerletListHybridTabulated(
            verletlistCG, True
        ))
    if lj_interaction is not None:
        system.addInteraction(lj_interaction, 'xyz-lj')
    if coulomb_interaction is not None:
        system.addInteraction(coulomb_interaction, 'xyz-coulomb')
    if tab_cg is not None:
        system.addInteraction(tab_cg, 'xyz-cg')

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
    coulomb_interaction = espressopp.interaction.VerletListHybridReactionFieldGeneralized(
        verletlistAT, False)
    coulomb_interaction = gromacs_topology.setCoulombInteractions(
        system, verletlistAT, args.coulomb_cutoff, input_conf.atomtypeparams,
        epsilon1=args.coulomb_epsilon1,
        epsilon2=args.coulomb_epsilon2, kappa=args.coulomb_kappa,
        interaction=coulomb_interaction)
    pair14_interactions = tools.setPairInteractions(
        system, input_conf, args.lj_cutoff)
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
    tools.setBondedInteractions(
        system, input_conf, force_static=True, only_at=True)
    tools.setAngleInteractions(
        system, input_conf, force_static=True, only_at=True)
    tools.setDihedralInteractions(
        system, input_conf, force_static=True, only_at=True)

    return verletlistAT, verletlistCG
