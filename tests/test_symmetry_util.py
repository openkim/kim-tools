#!/usr/bin/python

import math

import numpy as np
from ase.build import bulk
from ase.calculators.kim.kim import KIM
from ase.io import read
from elastic_excerpt import get_unique_components_and_reconstruct_matrix
from sympy import matrix2numpy, symbols

from kim_tools import (
    CENTERING_DIVISORS,
    change_of_basis_atoms,
    get_change_of_basis_matrix_to_conventional_cell_from_formal_bravais_lattice,
    get_crystal_structure_from_atoms,
    get_formal_bravais_lattice_from_space_group,
    get_space_group_number_from_prototype,
    get_symbolic_cell_from_formal_bravais_lattice,
)
from kim_tools.symmetry_util.core import (
    PeriodExtensionException,
    fit_voigt_tensor_to_cell_and_space_group,
    kstest_reduced_distances,
    reduce_and_avg,
)


def test_change_of_basis_atoms(
    atoms_conventional=bulk("SiC", "zincblende", 4.3596, cubic=True)
):
    calc = KIM("LennardJones612_UniversalShifted__MO_959249795837_003")
    atoms_conventional.calc = calc
    crystal_structure = get_crystal_structure_from_atoms(atoms_conventional)
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


def test_fit_voigt_tensor_to_cell_and_space_group():
    a, b, c, alpha, beta, gamma = symbols("a b c alpha beta gamma")
    # taken from A5B11CD8E_aP26_1_5a_11a_a_8a_a-001
    test_substitutions = [
        (a, 1.0),
        (b, 1.20466246551),
        (c, 1.81123604761),
        (alpha, math.radians(76.515)),
        (beta, math.radians(81.528)),
        (gamma, math.radians(71.392)),
    ]
    # Generate a random symmetric matrix
    c = np.random.random((6, 6))
    c = c + c.T
    for sgnum in range(1, 231):
        lattice = get_formal_bravais_lattice_from_space_group(sgnum)
        symbolic_cell = get_symbolic_cell_from_formal_bravais_lattice(lattice)
        cell = matrix2numpy(symbolic_cell.subs(test_substitutions), dtype=float)
        # This takes any matrix, picks out the unique constants based on the
        # algebraic diagrams, and returns a matrix conforming to the material symmetry

        # Test 1: fit_voigt_tensor_to_cell_and_space_group should not change
        # a matrix that is already symmetrized
        _, _, c_sym_alg = get_unique_components_and_reconstruct_matrix(c, sgnum)
        c_sym_alg_tens = fit_voigt_tensor_to_cell_and_space_group(
            c_sym_alg, cell, sgnum
        )
        assert np.allclose(c_sym_alg, c_sym_alg_tens)

        # Test 2: fit_voigt_tensor_to_cell_and_space_group should produce
        # a matrix that is already symmetrized
        c_sym_tens = fit_voigt_tensor_to_cell_and_space_group(c, cell, sgnum)
        if lattice == "aP":
            assert np.allclose(c, c_sym_tens)
        _, _, c_sym_tens_alg = get_unique_components_and_reconstruct_matrix(
            c_sym_tens, sgnum
        )
        assert np.allclose(c_sym_tens, c_sym_tens_alg)


if __name__ == "__main__":
    test_fit_voigt_tensor_to_cell_and_space_group()
