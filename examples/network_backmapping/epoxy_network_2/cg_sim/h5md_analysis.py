# Pierre de Buyl 2014
# This file is licensed under the modified BSD license

import espressopp

import numpy as np  # NOQA
import pyh5md


class DumpH5MD(object):
    def __init__(self, filename, system, integrator, edges, particle_ids,
                 save_vel=False, save_force=False, unfolded=False):
        self.system = system
        self.integrator = integrator

        self.unfolded = unfolded
        self.save_vel = save_vel
        self.save_force = save_force

        self.N = len(particle_ids)
        self.particle_ids = particle_ids
        self.h5file = pyh5md.H5MD_File(
            filename,
            'w',
            creator='espressopppp',
            creator_version=espressopp.Version().info(),
            author='XXX',
            author_email='XXX'
        )
        self.atoms_ds = self.h5file.particles_group('atoms')
        self.pos_ds = self.atoms_ds.trajectory('position', (self.N, 3), np.float64)
        self.image_ds = self.atoms_ds.trajectory('image', (self.N, 3), np.float64)
        # Saves particle types.
        self.atoms_ds.trajectory(
            'species',
            (self.N, ),
            np.int32,
            data=np.array([system.storage.getParticle(pid).type for pid in particle_ids]),
            time=False
        )
        self.atoms_ds.trajectory(
            'mass',
            (self.N, ),
            np.float64,
            data=np.array([system.storage.getParticle(pid).mass for pid in particle_ids]),
            time=False)
        if save_force:
            self.force_ds = self.atoms_ds.trajectory('force', (self.N, 3), np.float64)
        if save_vel:
            self.vel_ds = self.atoms_ds.trajectory('velocity', (self.N, 3), np.float64)
        self.h5file.box = self.atoms_ds.box(
            dimension=3,
            boundary=['periodic', 'periodic', 'periodic'],
            edges=edges)
        self.h5file.NPart = espressopp.analysis.NPart(system).compute()
        self.h5file.observable('particle_number', data=int(self.h5file.NPart), time=False)
        observables = [
            ('temperature', (), np.float64),
            ('kinetic_energy', (), np.float64),
            ('pressure', (), np.float64),
            ('pressure_tensor', (6,), np.float64),
            ('potential_energy', (), np.float64),
            ('total_energy', (), np.float64),
            ('resolution', (), np.float64)
        ]
        self.obs_dict = {
            o[0]: self.h5file.observable(*o) for o in observables
        }
        # Puts the interactions.
        for key in system.getInteractionLabels():
            self.obs_dict[key] = self.h5file.observable(key, (), np.float64)

    def dump(self):
        step = self.integrator.step
        time = self.integrator.step*self.integrator.dt
        data = []
        image = []
        for pid in self.particle_ids:
            p = self.system.storage.getParticle(pid)
            data.append(list(p.pos))
            image.append(list(p.imageBox))
        self.pos_ds.append(data, step, time)
        self.image_ds.append(image, step, time)
        
        if self.save_vel:
            self.vel_ds.append(
                [[x for x in self.system.storage.getParticle(pid).v] for pid in self.particle_ids],
                step,
                time
            )
        if self.save_force:
            self.force_ds.append(
                [[x for x in self.system.storage.getParticle(pid).f] for pid in self.article_ids],
                step,
                time
            )

    def flush(self):
        self.h5file.flush()

    def close(self):
        self.h5file.close()

    def analyse(self):
        step = self.integrator.step
        time = self.integrator.step*self.integrator.dt
        T = espressopp.analysis.Temperature(self.system).compute() / self.system.kb
        self.obs_dict['temperature'].append(T, step, time)
        P = espressopp.analysis.Pressure(self.system).compute()
        self.obs_dict['pressure'].append(P, step, time)
        Pij = espressopp.analysis.PressureTensor(self.system).compute()
        self.obs_dict['pressure_tensor'].append(Pij, step, time)

        Ek = (3.0/2.0) * T
        self.obs_dict['kinetic_energy'].append(Ek, step, time)

        # Gets the potential energy
        potential_energy = 0.0
        for interaction_label, interaction_id in self.system.getInteractionLabels().iteritems():
            energy_value = self.system.getInteraction(interaction_id).computeEnergy()
            self.obs_dict[interaction_label].append(energy_value, step, time)
            potential_energy += energy_value

        self.obs_dict['potential_energy'].append(potential_energy, step, time)
        self.obs_dict['total_energy'].append(Ek + potential_energy, step, time)


