import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.lj import LennardJones, cutoff_function, d_cutoff_function
from ase.neighborlist import NeighborList
from ase.stress import full_3x3_to_voigt_6_stress


class LennardJonesFailNoNeighbors(LennardJones):
    """
    A clone of ASE Lennard-Jones that raises an error
    if there are no neighbors
    """

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        sigma = self.parameters.sigma
        epsilon = self.parameters.epsilon
        rc = self.parameters.rc
        ro = self.parameters.ro
        smooth = self.parameters.smooth

        if self.nl is None or "numbers" in system_changes:
            self.nl = NeighborList(
                [rc / 2] * natoms, self_interaction=False, bothways=True
            )

        self.nl.update(self.atoms)

        positions = self.atoms.positions
        cell = self.atoms.cell

        # potential value at rc
        e0 = 4 * epsilon * ((sigma / rc) ** 12 - (sigma / rc) ** 6)

        energies = np.zeros(natoms)
        forces = np.zeros((natoms, 3))
        stresses = np.zeros((natoms, 3, 3))

        for ii in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(ii)

            # Testing raising an error for empty neighbor list
            assert len(neighbors) > 0

            cells = np.dot(offsets, cell)

            # pointing *towards* neighbours
            distance_vectors = positions[neighbors] + cells - positions[ii]

            r2 = (distance_vectors**2).sum(1)
            c6 = (sigma**2 / r2) ** 3
            c6[r2 > rc**2] = 0.0
            c12 = c6**2

            if smooth:
                cutoff_fn = cutoff_function(r2, rc**2, ro**2)
                d_cutoff_fn = d_cutoff_function(r2, rc**2, ro**2)

            pairwise_energies = 4 * epsilon * (c12 - c6)
            pairwise_forces = -24 * epsilon * (2 * c12 - c6) / r2  # du_ij

            if smooth:
                # order matters, otherwise the pairwise energy is already
                # modified
                pairwise_forces = (
                    cutoff_fn * pairwise_forces + 2 * d_cutoff_fn * pairwise_energies
                )
                pairwise_energies *= cutoff_fn
            else:
                pairwise_energies -= e0 * (c6 != 0.0)

            pairwise_forces = pairwise_forces[:, np.newaxis] * distance_vectors

            energies[ii] += 0.5 * pairwise_energies.sum()  # atomic energies
            forces[ii] += pairwise_forces.sum(axis=0)

            stresses[ii] += 0.5 * np.dot(
                pairwise_forces.T, distance_vectors
            )  # equivalent to outer product

        # no lattice, no stress
        if self.atoms.cell.rank == 3:
            stresses = full_3x3_to_voigt_6_stress(stresses)
            self.results["stress"] = stresses.sum(axis=0) / self.atoms.get_volume()
            self.results["stresses"] = stresses / self.atoms.get_volume()

        energy = energies.sum()
        self.results["energy"] = energy
        self.results["energies"] = energies

        self.results["free_energy"] = energy

        self.results["forces"] = forces
