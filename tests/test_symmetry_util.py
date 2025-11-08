#!/usr/bin/python

import math

import numpy as np
from ase.build import bulk
from ase.calculators.kim.kim import KIM
from ase.geometry import get_duplicate_atoms
from ase.io import read
from elastic_excerpt import algebraically_reconstruct_matrix
from sympy import matrix2numpy, symbols

from kim_tools import (
    CENTERING_DIVISORS,
    change_of_basis_atoms,
    get_atoms_from_crystal_structure,
    get_change_of_basis_matrix_to_conventional_cell_from_formal_bravais_lattice,
    get_crystal_structure_from_atoms,
    get_formal_bravais_lattice_from_space_group,
    get_primitive_genpos_ops,
    get_space_group_number_from_prototype,
    get_symbolic_cell_from_formal_bravais_lattice,
)

# Disorganized to test it here, but the operations are really redundant
# with testing transform_atoms
from kim_tools.aflow_util.core import get_atom_indices_for_each_wyckoff_orb
from kim_tools.symmetry_util.core import (
    FixProvidedSymmetry,
    PeriodExtensionException,
    fit_voigt_tensor_and_error_to_cell_and_space_group,
    fit_voigt_tensor_to_cell_and_space_group,
    reduce_and_avg,
    transform_atoms,
)


def test_change_of_basis_atoms(
    atoms_conventional=bulk("SiC", "zincblende", 4.3596, cubic=True)
):
    calc = KIM("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
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
        "structures/FeP_unstable.extxyz": True,
        "structures/FeP_stable.extxyz": False,
    }
    repeat = [10, 10, 10]
    for data_file in data_file_has_period_extension:
        has_period_extension = data_file_has_period_extension[data_file]
        atoms = read(data_file)
        try:
            reduce_and_avg(atoms, repeat)
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

    ########################################################
    # full coverage: add 3, 16, 142, 143, 152, 153, 168, 195
    # These tests are very slow, ~30sec each. So just test
    # a random triclinic matrix (to check that it doesn't
    # change under symmetrization), and a specific tetragonal
    # matrix with error.
    ########################################################

    for sgnum in (1,):
        lattice = get_formal_bravais_lattice_from_space_group(sgnum)
        symbolic_cell = get_symbolic_cell_from_formal_bravais_lattice(lattice)
        cell = matrix2numpy(symbolic_cell.subs(test_substitutions), dtype=float)
        # Throw some exponents in there
        cell += np.random.rand(3, 3) * 1e-16

        # Not checking the error here, so just pass c itself as a dummy error
        c_mat_symm_rot_sympy, _ = fit_voigt_tensor_and_error_to_cell_and_space_group(
            c, c, cell, sgnum
        )
        c_mat_symm_rot = fit_voigt_tensor_to_cell_and_space_group(c, cell, sgnum)
        if lattice == "aP":
            assert np.allclose(c, c_mat_symm_rot_sympy)
            assert np.allclose(c, c_mat_symm_rot)
        # This takes any matrix, picks out the unique constants based on the
        # algebraic diagrams, and returns a matrix conforming to the material symmetry
        c_mat_symm_alg = algebraically_reconstruct_matrix(c_mat_symm_rot_sympy, sgnum)
        assert np.allclose(c_mat_symm_rot_sympy, c_mat_symm_alg)
        assert np.allclose(c_mat_symm_rot, c_mat_symm_alg)

    # Use SG 75 to test a specific matrix including error
    sgnum = 75
    lattice = get_formal_bravais_lattice_from_space_group(sgnum)
    symbolic_cell = get_symbolic_cell_from_formal_bravais_lattice(lattice)
    cell = matrix2numpy(symbolic_cell.subs(test_substitutions), dtype=float)
    # Throw some exponents in there
    cell += np.random.rand(3, 3) * 1e-16

    c = np.ones((6, 6))
    c_err = np.ones((6, 6))

    # In Laue Class 4/m, [0,5]=-[1,5].
    # Here we have [0,5]=0, [1,5]=1, so
    # we should end up with [0,5]=-1/2,
    # [1,5]=1/2, and the resulting variance
    # be 1/4 of the sum of their variances
    c[0, 5] = 0
    c[5, 0] = 0
    c_err[1, 5] = 15 ** (1 / 2)
    c_err[5, 1] = 15 ** (1 / 2)

    ref_out = (
        [
            [1.0, 1.0, 1.0, 0.0, 0.0, -0.5],
            [1.0, 1.0, 1.0, 0.0, 0.0, 0.5],
            [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [-0.5, 0.5, 0.0, 0.0, 0.0, 1.0],
        ],
        [
            [0.5**0.5, 0.5**0.5, 0.5**0.5, 0.0, 0.0, 2.0],
            [0.5**0.5, 0.5**0.5, 0.5**0.5, 0.0, 0.0, 2.0],
            [0.5**0.5, 0.5**0.5, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5**0.5, 0.5**0.5, 0.0],
            [0.0, 0.0, 0.0, 0.5**0.5, 0.5**0.5, 0.0],
            [2.0, 2.0, 0.0, 0.0, 0.0, 1.0],
        ],
    )
    ref_out_symm = (
        [
            [1.0, 1.0, 1.0, 0.0, 0.0, -0.5],
            [1.0, 1.0, 1.0, 0.0, 0.0, 0.5],
            [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [-0.5, 0.5, 0.0, 0.0, 0.0, 1.0],
        ],
        [
            [0.5**0.5, 1, 0.5**0.5, 0.0, 0.0, 2.0],
            [1, 0.5**0.5, 0.5**0.5, 0.0, 0.0, 2.0],
            [0.5**0.5, 0.5**0.5, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5**0.5, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5**0.5, 0.0],
            [2.0, 2.0, 0.0, 0.0, 0.0, 1.0],
        ],
    )
    assert np.allclose(
        fit_voigt_tensor_and_error_to_cell_and_space_group(c, c_err, cell, sgnum),
        ref_out,
    )
    assert np.allclose(
        fit_voigt_tensor_and_error_to_cell_and_space_group(c, c_err, cell, sgnum, True),
        ref_out_symm,
    )
    assert np.allclose(
        fit_voigt_tensor_to_cell_and_space_group(c, cell, sgnum),
        ref_out[0],
    )


def test_FixProvidedSymmetry():
    test_cases = [
        {
            "symm": 1,
            "deform_allowed": True,
        },
        {
            "symm": get_primitive_genpos_ops(1),
            "deform_allowed": True,
        },
        {
            "symm": 221,
            "deform_allowed": False,
        },
        {
            "symm": get_primitive_genpos_ops(221),
            "deform_allowed": False,
        },
    ]
    for case in test_cases:
        atoms = bulk("CsCl", "cesiumchloride", 1.0)
        constraint = FixProvidedSymmetry(atoms, case["symm"])
        atoms.set_constraint(constraint)
        atoms.set_positions(
            [
                [
                    0,
                    0,
                    0,
                ],
                [0.6, 0.5, 0.5],
            ]
        )
        deform_allowed = case["deform_allowed"]
        assert np.allclose(atoms[1].position, [0.6, 0.5, 0.5]) == deform_allowed
        new_cell = [[1.1, 0, 0], [0, 1, 0], [0, 0, 1]]
        atoms.set_cell(new_cell, scale_atoms=True)
        assert np.allclose(atoms.cell, new_cell) == deform_allowed


def test_transform_atoms_and_symmetry_and_wyckoff_indices(input_crystal_structures):
    for material in input_crystal_structures:
        atoms = get_atoms_from_crystal_structure(material)
        prototype_label = material["prototype-label"]["source-value"]
        wyckoff_atom_indices = get_atom_indices_for_each_wyckoff_orb(prototype_label)
        sgnum = get_space_group_number_from_prototype(prototype_label)
        genpos_ops = get_primitive_genpos_ops(sgnum)
        # Check that transforming by a symmetry op always gets an identical
        # atoms object
        for op in genpos_ops:
            atoms_transformed = transform_atoms(atoms, op)
            atoms_combined = atoms_transformed + atoms
            dupes = get_duplicate_atoms(atoms_combined, delete=True)
            assert len(atoms_combined) == len(atoms)
            # Check that any duplicates correspond to the same Wyckoff orbit
            for dupe in dupes:
                for position in wyckoff_atom_indices:
                    if dupe[0] in position["indices"]:
                        assert dupe[1] - len(atoms) in position["indices"]


if __name__ == "__main__":
    test_fit_voigt_tensor_to_cell_and_space_group()
