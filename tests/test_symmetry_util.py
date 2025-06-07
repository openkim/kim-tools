#!/usr/bin/python

import numpy as np
from ase.build import bulk
from ase.calculators.kim.kim import KIM
from ase.io import read

from kim_tools import (
    CENTERING_DIVISORS,
    change_of_basis_atoms,
    get_change_of_basis_matrix_to_conventional_cell_from_formal_bravais_lattice,
    get_crystal_structure_from_atoms,
    get_formal_bravais_lattice_from_space_group,
    get_space_group_number_from_prototype,
)
from kim_tools.symmetry_util.core import (
    PeriodExtensionException,
    kstest_reduced_distances,
    reduce_and_avg,
)


def test_change_of_basis_atoms(
    atoms_conventional=bulk("SiC", "zincblende", 4.3596, cubic=True)
):
    calc = KIM("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
    atoms_conventional.calc = calc
    crystal_structure = get_crystal_structure_from_atoms(
        atoms_conventional, get_short_name=False
    )
    prototype_label = crystal_structure["prototype-label"]["source-value"]
    sgnum = get_space_group_number_from_prototype(prototype_label)
    formal_bravais_lattice = get_formal_bravais_lattice_from_space_group(sgnum)
    primitive_to_conventional_change_of_basis = (
        get_change_of_basis_matrix_to_conventional_cell_from_formal_bravais_lattice(
            formal_bravais_lattice
        )
    )
    conventional_to_primitive_change_of_basis = np.linalg.inv(
        primitive_to_conventional_change_of_basis
    )
    centering = formal_bravais_lattice[1]
    multiplier = np.linalg.det(primitive_to_conventional_change_of_basis)
    assert np.isclose(multiplier, CENTERING_DIVISORS[centering])
    conventional_energy = atoms_conventional.get_potential_energy()
    atoms_primitive = change_of_basis_atoms(
        atoms_conventional, conventional_to_primitive_change_of_basis
    )
    atoms_primitive.calc = calc
    primitive_energy = atoms_primitive.get_potential_energy()
    assert np.isclose(primitive_energy * multiplier, conventional_energy)
    atoms_conventional_rebuilt = change_of_basis_atoms(
        atoms_primitive, primitive_to_conventional_change_of_basis
    )
    atoms_conventional_rebuilt.calc = calc
    conventional_rebuilt_energy = atoms_conventional_rebuilt.get_potential_energy()
    assert np.isclose(conventional_energy, conventional_rebuilt_energy)


def test_test_reduced_distances():
    data_file_has_period_extension = {
        "structures/FeP_period_extended_phase_transition.data": True,
        "structures/FeP_stable.data": False,
    }
    repeat = [11, 11, 11]
    for data_file in data_file_has_period_extension:
        has_period_extension = data_file_has_period_extension[data_file]
        atoms = read(data_file, format="lammps-data")
        _, reduced_distances = reduce_and_avg(atoms, repeat)
        try:
            kstest_reduced_distances(reduced_distances)
            assert not has_period_extension
        except PeriodExtensionException:
            assert has_period_extension


if __name__ == "__main__":
    test_change_of_basis_atoms()
